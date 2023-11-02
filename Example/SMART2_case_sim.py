import os
import PEEP_Ppleateau_detect_lib as PPD
import OptimalPEEP_lib as OP

directory_pressure = 'D:/LRZ Sync+Share/austausch_wei_teng_kao (Carolin Geitner)/data_simulation/SMART_E2'
filename_pressure = 'test_SMART_E2_lung_new_V15000_awac_nlin_data_Paw.csv'
csv_file_path_pressrue = os.path.join(directory_pressure, filename_pressure)

directory_flow = 'D:/LRZ Sync+Share/austausch_wei_teng_kao (Carolin Geitner)/data_simulation/SMART_E2'
filename_flow = 'test_SMART_E2_lung_new_V15000_awac_nlin_data_Q.csv'
csv_file_path_flow = os.path.join(directory_flow, filename_flow)

#Import the pressure and flow data
start = 816  
end = 995
result_P, index_P = PPD.store_csv_column(csv_file_path_pressrue, detect_label= 'avg(pressure)', initial_time = 0, start_time = start, end_time = end, plot_=False)
result_F, index_F = PPD.store_csv_column(csv_file_path_flow, detect_label= 'avg(flow_in)', initial_time = 0, start_time = start, end_time = end, plot_=False)

#Find Peak value
patient2 = PPD.ExtractPVDandPressure(result_P,result_F,frquency_=25)
pos,neg = patient2.find_peak( rm_unnecessary_pressure= True, time_over = 4,plot_= True, plot_all_= False)

#Only use for generating the essential files for holo-shed algorithm
patient2.createXMLfile(
    "D:/LRZ Sync+Share/austausch_wei_teng_kao (Carolin Geitner)/data_simulation/SMART_E2/test_SMART_E2_lung_new_V15000_awac_nlin_PEEPtrial/test_SMART_E2_lung_new_V15000_awac_nlin_PEEPtrial/test_SMART_E2_lung_new_V15000_awac_nlin-files",
    "D:\LRZ Sync+Share/austausch_wei_teng_kao (Carolin Geitner)/HIWI_WeiTengPythonCode/EIT_OptimalPEEP",
    "config_file_sim_template.json",
    "D:/LRZ Sync+Share/austausch_wei_teng_kao (Carolin Geitner)/data_simulation/SMART_E2/test_SMART_E2_lung_new_V15000_awac_nlin_PEEPtrial/test_SMART_E2_lung_new_V15000_awac_nlin_PEEPtrial/output_folder", 
    "output_SMART2.pvd"
    ,shift=0
    ,createPVD = False
    ,transfer2VTU= False    # if true, combine the  PVD files from all each post-preocssor in one VTU file
    ,createConfig_sim= True   
    ,createConfig_meas = False
    ,limit_sim = [20360, 24860]
    ,limit_meas = [0,-1])

# Final Step: finding the optimal PEEP
P_plateau = PPD.calculate_decremental_pressure(pos , NumValueCounted = -1, abandon_last=False)
PEEP = PPD.calculate_decremental_pressure(neg , NumValueCounted = -1, abandon_last=False)

dir = "D:/LRZ Sync+Share/austausch_wei_teng_kao (Carolin Geitner)/HIWI_WeiTengPythonCode/EIT_OptimalPEEP"
file_path_dict = {
    'level_1': os.path.join(dir,"output_EIT_1.txt"),
    'level_2': os.path.join(dir,"output_EIT_2.txt"),
    'level_3': os.path.join(dir,"output_EIT_3.txt"),
    'level_4': os.path.join(dir,"output_EIT_4.txt"),
    'level_5': os.path.join(dir,"output_EIT_5.txt"),
    'level_6': os.path.join(dir,"output_EIT_6.txt"),
}

OP.findOptimal(file_path_dict ,P_plateau, PEEP,'avg',-1,True )



