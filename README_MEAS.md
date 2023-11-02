# Project - Finding Optimal PEEP from measurement data
Author: Wei-Teng Kao, ge32gak@mytum.de

## Required library
- Open source:
1. vtk
2. numpy
3. pandas
4. scipy

Use the package manager [pip](https://pip.pypa.io/en/stable/) to install open source library.

```bash
pip install vtk
```
- self-written library
1. PEEP_Ppleateau_detect_lib
2. OptimalPEEP_lib


## Flow Chart
### 1. Importing .csv file. 
In order to find the PEEP and P-Plateau for decremental PEEP trial, pressure and flow data are required   

```python
import PEEP_Ppleateau_detect_lib as PPD

# pressure data
directory_pressure = 'D:/LRZ Sync+Share/austausch_wei_teng_kao (Carolin Geitner)/data/SMART_E2/SMART-E2 Studienprotokoll Tag1/csv/realtimeData/50Hz'
filename_pressure = 'realtimeValues_20210707094239_time.csv'
csv_file_path_pressrue = os.path.join(directory_pressure, filename_pressure)

# flow data
directory_flow = 'D:/LRZ Sync+Share/austausch_wei_teng_kao (Carolin Geitner)/data/SMART_E2/SMART-E2 Studienprotokoll Tag1/csv/realtimeData/50Hz'
filename_flow = 'realtimeValues_20210707094239_time.csv'
csv_file_path_flow = os.path.join(directory_flow, filename_flow)

# giving the bound of start_time and end_time
# Attension: in order to find the correct Peak, please select the start_time and end_time carefully, detailed explanation with the image will be below 
#Import the pressure and flow data
start = 1000
end = 1309

result_P, index_P = PPD.store_csv_column(csv_file_path_pressrue, detect_label='Pressure [mbar]', initial_time = 0, start_time = start, end_time = end,frquency = 50, plot_=False)
result_F, index_F = PPD.store_csv_column(csv_file_path_flow, detect_label='Flow [l/min]' , initial_time = 0, start_time = start, end_time = end,frquency = 50, plot_=False)

```

### 2. Finding the Peaks
- Attension! the recorded frequency of measurement data is 50Hz instead of 25Hz.
```python
# import the pressure data and flow data, and initialize the ExtractPVDandPressure object. 
patient2 = PPD.ExtractPVDandPressure(result_P,result_F,frquency_ = 50)

# Finding the peaks
# rm_unnecessary_pressure = True: it will remove the peaks of pressure between two different pressure levels.
# time_over = 3.5: if rm_unnecessary_pressure = True, it will remove the front peak if the interval of two peaks is smaller than 4
# plot_: plot only the non-scaled pressure w.r.t index
# plot_all_: plot all the information, including flow pressure. Attention: the value here are all normalized. 
pos,neg = patient2.find_peak( rm_unnecessary_pressure= True, time_over = 3.5,plot_= True, plot_all_= False)
# pos and neg represented the upper peals and lower peaks respectively. 
```

![Semantic description of image](/Pic_readme_used/Figure_1.png "Peaks dectection")

1. In order to find the correct peaks for simulation, the first peak of each pressure level **has to** be loacted at the PEEP,not the P-Plateau, meaning the started time frame should be choosed wisely between the P-plateau and PEEP. For Tthe last peak, it should be loacted at the P-plateau, so the ended time frame wll be between the PEEP and P_plateau. Here in the example, ```start_time = 1000, end_time = 1309```.
2. The peaks denoted by star symbol are the peaks will be used for further analysis. e.g. In the next step, Create the necessary files for the holo-shed algorithm, the index of those peaks are selected for generating the ```.json```files. 

### 3. Create the necessary files for the holo-shed algorithm
1. First, it will generate the .pvd file. you can import this file into PARAVIEW to visualize the results and check the position of the detected peak is correct or not, if you face a synchronized problem when comparing imported data and detected Peaks, you can also use this way to calculate the time-shifting manually.
2. if ```transferVTU = True```: it will also generate the ```VTU``` file. These are necessary for running holo-shed algorithm. **ATTENTION: please make sure all the index in the range has their corresponding ```pVTU``` files, otherwise the error will occur.**
3. if ```createConfig_meas = True```: it will generate the ```config_file.json```. These are necessary for running holo-shed algorithm.
4. **This function is only used when you need to generate the essential files for holo-shed algorithm, if you already have them, then you can ignore this function and comment this function, there is no influence on finding the optimal PEEP** 
5. limit_meas and limit_sim determine the range of the export data. Notice that for simulation case, the range is determined by frame index and for measurement case, it is determined by time (sec).

```python.
#Only use for generating the essential files for holo-shed algorithm
patient2.createXMLfile(
    "D:/LRZ Sync+Share/austausch_wei_teng_kao (Carolin Geitner)/data_simulation/SMART_E2/test_SMART_E2_lung_new_V15000_awac_nlin_PEEPtrial/test_SMART_E2_lung_new_V15000_awac_nlin_PEEPtrial/test_SMART_E2_lung_new_V15000_awac_nlin-files",
    "D:/LRZ Sync+Share/austausch_wei_teng_kao (Carolin Geitner)/data_simulation/SMART_E2/test_SMART_E2_lung_new_V15000_awac_nlin_PEEPtrial/test_SMART_E2_lung_new_V15000_awac_nlin_PEEPtrial/output_folder", 
    "output_SMART2.pvd"
    ,shift=0
    ,createPVD = False
    ,transfer2VTU= False    # if true, combine the  PVD files from all each post-preocssor in one VTU file
    ,createConfig_sim= False   
    ,createConfig_meas = True
    ,limit_sim = [20360, 24860]
    ,limit_meas = [0,-1])
```

### 4. Generate the EIT data from the holo-shed algorithm
1. **ATTENTION, This is done using the holo-shed algorithm, please check [holoshed](for the instruction**) for the instruction**
2. **In order to extract the simulated EIT value for further tasks, please add these lines in the ```main.py```of holo-shed algorithm, It can help to store the EIT value in .txt files**
3. Version of holoshed-algorithm: ```1a317a4b435efc511bd41d64a7f7ff25902717a8```,from branch: ```enable-reconstruction-of-measurements-only```
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
        EIT_meas_array = np.array(EIT_meas)
        np.savetxt('output_EIT_meas_9.txt',EIT_meas_array.reshape(len(EIT_meas_array),-1),delimiter=',')
        holoshed.plot_global.plot_detached_electrodes(config, interactive)
        #############################################################################################
        #############################################################################################
        #############################################################################################

    holoshed.plot_global.plot_global(config, interactive)
    print("Execution of Holo-Shed successfull.")
    return 0


if __name__ == "__main__":
    main()
```
After this, you will obtain several results w.r.t each pressure level, They will be stored as ```.txt```format. It should be a 2D-array where each row represents one flatten EIT image.

### 5. Compute the Optimal PEEP
1. Back to our own code, in order to compute the Optimal PEEP, you will need the PEEP and P-plateau value. They can be extracted by calling ```calculate_decremental_pressure()```
```python
# Final Step: finding the optimal PEEP
P_plateau = PPD.calculate_decremental_pressure(pos ,p_thershold= 1, NumValueCounted = -1, abandon_last=False)
PEEP = PPD.calculate_decremental_pressure(neg,p_thershold= 1, NumValueCounted = -1, abandon_last=False)
```
2. After that, create the dictionary for the EIT array that is computed by the holo-shed algorithm,
it will be used when calling  
```OP.findOptimal()```function
```python
#import the library for finding Optimal PEEP
import OptimalPEEP_lib as OP

#Formulate the dict of EIT array
dir = "D:/LRZ Sync+Share/austausch_wei_teng_kao (Carolin Geitner)/HIWI_WeiTengPythonCode/EIT_OptimalPEEP"
file_path_dict = {
    'level_0': os.path.join(dir,"output_EIT_meas_0.txt"),
    'level_1': os.path.join(dir,"output_EIT_meas_1.txt"),
    'level_2': os.path.join(dir,"output_EIT_meas_2.txt"),
    'level_3': os.path.join(dir,"output_EIT_meas_3.txt"),
    'level_4': os.path.join(dir,"output_EIT_meas_4.txt"),
    'level_5': os.path.join(dir,"output_EIT_meas_5.txt"),
    'level_6': os.path.join(dir,"output_EIT_meas_6.txt"),
    'level_7': os.path.join(dir,"output_EIT_meas_7.txt"),
    'level_8': os.path.join(dir,"output_EIT_meas_8.txt"),
    'level_9': os.path.join(dir,"output_EIT_meas_9.txt"),
}

#calling the OP.findOptimal() function to compute the final result
OP.findOptimal(file_path_dict ,P_plateau*100, PEEP*100,'avg',-1,True )

```

## Example script
Here the whole script is presented, which is refer to SMART2 case.

First, please prepare the required files for the algorithm
1. For finding Peaks
- pressure data (```.csv``` file)
- flow data (```.csv``` file)
2. For ```CreateXML()``` function: Create ```.xml```, ```.pvtu```, ```.vtu```,```.json``` files for holo-shed algorithm
- Corresponding ```.pvtu``` data
- Corresponding ```.vtu``` data
3. For the holo-shed algorithm, please check [holo-shed](https://gitlab.com/xenotaph3/holoshed) for further information
- Corresponding ```.eit``` data, please input to ```01_input``` folder
- Config file (```.json``` file), generated by ```CreateXML()``` function, please input to ```00_config``` folder
- Correspoding measured data (50Hz) (```.csv```file)
- Correspoding mesh file
4. For finding optimal PEEP
- output EIT data (```.csv``` file)of each pressure level, generated by holo-shed algorithm, please make sure you modify the ```main.py``` file as step 4.
- pressure data (```.csv``` file), same as the one for finding Peaks
- flow data (```.csv``` file), same as the one for finding Peaks


```python
import os
import PEEP_Ppleateau_detect_lib as PPD
import OptimalPEEP_lib as OP


#'D:/LRZ Sync+Share/austausch_wei_teng_kao (Carolin Geitner)/data_simulation/SMART_E2'
#'test_SMART_E2_lung_new_V15000_awac_nlin_data_Paw.csv'
directory_pressure = 'D:/LRZ Sync+Share/austausch_wei_teng_kao (Carolin Geitner)/data/SMART_E2/SMART-E2 Studienprotokoll Tag1/csv/realtimeData/50Hz'
filename_pressure = 'realtimeValues_20210707094239_time.csv'
csv_file_path_pressrue = os.path.join(directory_pressure, filename_pressure)


directory_flow = 'D:/LRZ Sync+Share/austausch_wei_teng_kao (Carolin Geitner)/data/SMART_E2/SMART-E2 Studienprotokoll Tag1/csv/realtimeData/50Hz'
filename_flow = 'realtimeValues_20210707094239_time.csv'
csv_file_path_flow = os.path.join(directory_flow, filename_flow)

#Import the pressure and flow data
start = 1000
end = 1309

result_P, index_P = PPD.store_csv_column(csv_file_path_pressrue, detect_label='Pressure [mbar]', initial_time = 0, start_time = start, end_time = end,frquency = 50, plot_=False)
result_F, index_F = PPD.store_csv_column(csv_file_path_flow, detect_label='Flow [l/min]' , initial_time = 0, start_time = start, end_time = end,frquency = 50, plot_=False)

#Find Peak value
patient2 = PPD.ExtractPVDandPressure(result_P,result_F,50)
pos,neg = patient2.find_peak( rm_unnecessary_pressure= True, time_over = 3.5,plot_= True, plot_all_= False)


#Only use for generating the essential files for holo-shed algorithm
patient2.createXMLfile(
    "D:/LRZ Sync+Share/austausch_wei_teng_kao (Carolin Geitner)/data_simulation/SMART_E2/test_SMART_E2_lung_new_V15000_awac_nlin_PEEPtrial/test_SMART_E2_lung_new_V15000_awac_nlin_PEEPtrial/test_SMART_E2_lung_new_V15000_awac_nlin-files",
    "D:/LRZ Sync+Share/austausch_wei_teng_kao (Carolin Geitner)/data_simulation/SMART_E2/test_SMART_E2_lung_new_V15000_awac_nlin_PEEPtrial/test_SMART_E2_lung_new_V15000_awac_nlin_PEEPtrial/output_folder", 
    "output_SMART2.pvd"
    ,shift=0
    ,createPVD = False
    ,transfer2VTU= False    # if true, combine the  PVD files from all each post-preocssor in one VTU file
    ,createConfig_sim= False   
    ,createConfig_meas = True
    ,limit_sim = [20360, 24860]
    ,limit_meas = [0,-1])

# Final Step: finding the optimal PEEP
P_plateau = PPD.calculate_decremental_pressure(pos ,p_thershold= 1, NumValueCounted = -1, abandon_last=False)
PEEP = PPD.calculate_decremental_pressure(neg,p_thershold= 1, NumValueCounted = -1, abandon_last=False)

dir = "D:/LRZ Sync+Share/austausch_wei_teng_kao (Carolin Geitner)/HIWI_WeiTengPythonCode/EIT_OptimalPEEP"

file_path_dict = {
    'level_0': os.path.join(dir,"output_EIT_meas_0.txt"),
    'level_1': os.path.join(dir,"output_EIT_meas_1.txt"),
    'level_2': os.path.join(dir,"output_EIT_meas_2.txt"),
    'level_3': os.path.join(dir,"output_EIT_meas_3.txt"),
    'level_4': os.path.join(dir,"output_EIT_meas_4.txt"),
    'level_5': os.path.join(dir,"output_EIT_meas_5.txt"),
    'level_6': os.path.join(dir,"output_EIT_meas_6.txt"),
    'level_7': os.path.join(dir,"output_EIT_meas_7.txt"),
    'level_8': os.path.join(dir,"output_EIT_meas_8.txt"),
    'level_9': os.path.join(dir,"output_EIT_meas_9.txt"),
}

OP.findOptimal(file_path_dict ,P_plateau*100, PEEP*100,'avg',-1,True )
```

## Reference:
Costa EL, Borges JB, Melo A, Suarez-Sipmann F, Toufen C Jr, Bohm SH, Amato MB. Bedside estimation of recruitable alveolar collapse and hyperdistension by electrical impedance tomography. Intensive Care Med. 2009 Jun;35(6):1132-7. doi: 10.1007/s00134-009-1447-y. Epub 2009 Mar 3. PMID: 19255741.
