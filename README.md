# HIWI - Optimal PEEP finding from EIT

**HIWI-Optimal PEEP finding from EIT** is the python code for the HIWI. The main task is to evaluate the recruitable alveolar
collapse and hyperdistension by electrical impedance tomography ([Costa, 2009](https://pubmed.ncbi.nlm.nih.gov/19255741/)) and find the optimal PEEP for the patient.
To generate the required EIT data for the evaluation, user have to run the [Holoshed alogorithm](https://gitlab.com/xenotaph3/holoshed) then extract the data manually and store as the ```.txt``` file.

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
- Self-written library
  1. PEEP_Ppleateau_detect_lib
  2. OptimalPEEP_lib

## Generate ```EIT.txt``` data from HoloShed algorithm
 - [EIT ReadMe](/doc/GenerateEIT_holoshed.md)

## Example
 - [SMART2 CASE](/Example/Example_for_meas.ipynb)
  
## Reference
 Costa EL, Borges JB, Melo A, Suarez-Sipmann F, Toufen C Jr, Bohm SH, Amato MB. Bedside estimation of recruitable alveolar collapse and hyperdistension by electrical impedance tomography. Intensive Care Med. 2009 Jun;35(6):1132-7. doi: 10.1007/s00134-009-1447-y. Epub 2009 Mar 3. PMID: 19255741.
  










