## Generate the EIT data from the holo-shed algorithm
1. ATTENTION, This is done using the holoshed algorithm, please check [holoshed](https://gitlab.com/xenotaph3/holoshed) for the instruction
   
2. **In order to extract the simulated EIT value for further tasks, please add following lines in the ```main.py```of holo-shed algorithm, It will store the EIT values of all the images in current pressure level as flatten matrix in ```.txt``` files**
   
3. Version of holoshed-algorithm: ```1a317a4b435efc511bd41d64a7f7ff25902717a8```, from branch: ```enable-reconstruction-of-measurements-only```
   
```python
"""Main module of holoshed."""
import os

import click

import holoshed.plot_eit
import holoshed.plot_fe
import holoshed.plot_global
import holoshed.utils
from holoshed.initialize_holoshed import (
    initialize_holoshed as initialize_holoshed,
)
from holoshed.lung_fusion import lungFusion
from holoshed.makeGREIT import makeGREIT
from holoshed.setup_forward_model import setup_forward_model
import numpy as np


@click.command()
@click.argument("config_file", required=True, type=click.Path(exists=True))
@click.option(
    "--interactive",
    is_flag=True,
    help="Interactive mode. Show plots for debugging purposes.",
)
@click.option(
    "--write-voltage",
    is_flag=True,
    help="Write simulated voltage output for every timestep into ensight format.",
)
@click.option(
    "--write-conductivity",
    is_flag=True,
    help="Defines whether to write the conductivity for debugging purpose.",
)
def main(config_file, interactive, write_voltage, write_conductivity):
    """Main function of holoshed.
    Args:
        interactive         (bool):     Run holoshed in interactive mode, i.e.
                                        plots are shown during runtime for
                                        debugging purposes.
        write_voltage       (bool):     Write simulated electrode voltages
                                        to ensight output.
        write_conductivity  (bool):     Defines whether to write the conductivity
                                        for debugging purpose
    """
    config = initialize_holoshed(config_file)

    fmdl = setup_forward_model(config)

    if not config["run_meas_data_only"]:
        config = lungFusion(
            config=config,
            fmdl=fmdl,
            write_voltage=write_voltage,
            write_conductivity=write_conductivity,
        )
    holoshed.plot_fe.plot_fe_surf(fmdl, config, interactive)

    sim_eit_img, meas_eit_img, config = makeGREIT(config, fmdl)

    if not config["run_meas_data_only"]:
        output_dir = os.path.join(
            config["working_path"], "04_simulated_reconstruction_results"
        )
        holoshed.plot_eit.plot_eit(
            sim_eit_img, config, output_dir, config["sim_name"], interactive
        )

    if not config["run_sim_data_only"]:
        output_dir = os.path.join(
            config["working_path"], "05_measured_reconstruction_results"
        )
        EIT_meas = holoshed.plot_eit.plot_eit(
            meas_eit_img, config, output_dir, "meas", interactive
        )
        #############################################################################################
        ########################Please add the lines below to store the data#########################
        #############################################################################################

        output_filename = 'output_EIT_meas_0.txt'

        EIT_meas_array = np.array(EIT_meas)
        np.savetxt(output_filename,EIT_meas_array.reshape(len(EIT_meas_array),-1),delimiter=',')
        
        #############################################################################################
        #############################################################################################
        #############################################################################################

        holoshed.plot_global.plot_detached_electrodes(config, interactive)
    holoshed.plot_global.plot_global(config, interactive)
    print("Execution of Holo-Shed successfull.")
    return 0


if __name__ == "__main__":
    main()
```
After this, user **have to** run it for each pressure level (remember to modify the output_filename every time, otherwise the ```.txt``` file will be overwrite). It should be a 2D-array where each row represents one flatten EIT image.

## Flow chart

1. Modify the ```main.py``` as above.
2. Run holoshed in terminal
```python 
python main.py /home/v196mp6/EIT_measure_only/holoshed/00_config/config_file_meas_0.json
```
3. Check the result plot in ```06_global_postprocessing_results``` folder. The result should look like the following, which the peak of pressure be selected correctly and the peak of voltage located at the same time/index of peak of pressure
   
![reuslt](/demo_imag/demo.png)

4 Modify the ```main.py``` again, change the output file name of ```.txt``` file as  ```output_EIT_meas_1.txt``` for next pressure level.
   
1. Run holoshed in terminal again for next pressure level
```python 
python main.py /home/v196mp6/EIT_measure_only/holoshed/00_config/config_file_meas_1.json
```
1. Repeat step 3 to step 5 again and again until finish the simulation for all the pressure level. In this case (SMART2), user have to run it 10 times.

