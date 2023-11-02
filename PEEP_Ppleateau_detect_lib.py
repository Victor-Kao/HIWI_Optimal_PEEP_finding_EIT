import csv
import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
from scipy.signal import find_peaks, peak_widths
import xml.etree.ElementTree as ET
from scipy import signal
import vtk
import json


def store_csv_column(csv_file, detect_label = " ", initial_time = 0,  start_time = 0.0, end_time = -1.0, frquency = 25 , plot_ = False):
    '''
    Input:
        csv_file: ```string```, input csv file
        initial_time: ```float```, started time from the csv files, not the started time of array
        start_time: ```float```,  when the array start
        end_time: ```float```,  when the array end, default value = -1, meaning all the value will be stored
        frequency: ```int```, the frquency of the importing data
        plot_ : ```bool```,  plot the info for not, where x indicate time and y indicate avg pressure.
    Output:
        reuslt: ```np.array```,  output array, where store the time and avg pressure
        index: ```np.array```
,
    '''
    if detect_label == " ":
        result = np.array([])
        index = np.array([])
    else:
        df = pd.read_csv(csv_file, usecols=['Time',detect_label])
        start_index = int(( start_time - initial_time )*int(frquency))
        end_index = int(( end_time - initial_time )*int(frquency))

        if end_index  == -1:
            df = df.iloc[start_index:len(df),:]
            index =  np.arange(start_index, len(df), 1)
        else:
            df = df.iloc[start_index:end_index,:]
            index =  np.arange(start_index, end_index, 1)
        result = np.array(df)
    
    
    if plot_:
        plt.plot(result[:,0] , result[:,1])
        plt.xlabel('Time')
        plt.ylabel('avg(pressure)')
        plt.title('Plot of CSV Column')
        plt.grid(True)
        plt.show()

    return result, index

def create_config_HoloShed_meas(time_list, template_dir,template_name, output_dir,output_file_name):
    
    realtime_start = time_list[0]
    realtime_end = time_list[-1]
    realtime_reference = realtime_start
    realtime_reconst_times = time_list[:-2]  #abandon the last PEEP and P-plateau value.

    output_file_ = os.path.join(output_dir, output_file_name)
    # Create a dictionary with the desired data

    template_path_ = os.path.join(template_dir,template_name)
    # Open the file in read mode
    with open(template_path_, 'r') as file:
        # Read the entire file
        file_contents = file.read()
        data_json = json.loads(file_contents)

        if 'realtime_start' in data_json:
            data_json['realtime_start'] = realtime_start
        else:
            print("ERROR! missing key: realtime_start, please make sure the key is included in the template file (meas)")
        
        if 'realtime_end' in data_json:
            data_json['realtime_end'] = realtime_end
        else:
            print("ERROR! missing key: realtime_end, please make sure the key is included in the template file (meas)")

        if 'realtime_reference' in data_json:
            data_json['realtime_reference'] = realtime_reference
        else:
            print("ERROR! missing key: realtime_reference, please make sure the key is included in the template file (meas)")

        if 'realtime_reconst_times' in data_json:
            data_json['realtime_reconst_times'] = realtime_reconst_times.tolist()
        else:
            print("ERROR! missing key: realtime_reconst_times, please make sure the key is included in the template file (meas)")

    with open(output_file_, 'w') as json_file:
        json.dump(data_json, json_file, indent=1)


def create_config_HoloShed_sim(start_timestep_, end_timestep_, template_dir ,template_name, output_dir,output_file_name):
    '''
    - the data created is for holo-shed algorithm
    - after generate the json file, please put it into 00_config folders
    - ATTENSION: this function only setting the start_timestep and end_timestep for simulation, please check holo-shed algorithm for other default parameters
    '''
    output_file_ = os.path.join(output_dir, output_file_name)

    #template_dir = "D:\LRZ Sync+Share/austausch_wei_teng_kao (Carolin Geitner)/HIWI_WeiTengPythonCode/EIT_OptimalPEEP"
    #template_name = "config_file_sim_template.json"

    template_path_ = os.path.join(template_dir,template_name)
    # Open the file in read mode
    with open(template_path_, 'r') as file:
        # Read the entire file
        file_contents = file.read()
        data_json = json.loads(file_contents)

        if 'start_timestep' in data_json:
            data_json['start_timestep'] = start_timestep_
        else:
            print("ERROR! missing key: realtime_start, please make sure the key is included in the template file (meas)")
        
        if 'end_timestep' in data_json:
            data_json['end_timestep'] = end_timestep_
        else:
            print("ERROR! missing key: realtime_end, please make sure the key is included in the template file (meas)")

        if 'reference_timestep' in data_json:
            data_json['reference_timestep'] = start_timestep_
        else:
            print("ERROR! missing key: realtime_reference, please make sure the key is included in the template file (meas)")

    # Write the dictionary to a JSON file
    with open(output_file_, 'w') as json_file:
        json.dump(data_json, json_file, indent=1)


class ExtractPVDandPressure:
    def __init__(self, input_array_pressure,input_array_flow, frquency_ ):
        self.time = input_array_pressure[:,0]
        self.pressure = input_array_pressure[:,1]
        self.pressure_index = input_array_pressure[:,0]
        self.flow = input_array_flow[:,1]
        self.frequency = frquency_

    def generateVTUdata (self, inputfile, output_dir, filename):
        '''
        - Generate the ```.vtu``` files by combining fractional```.vtu``` files generated from post-processor.
        '''
        reader = vtk.vtkXMLPUnstructuredGridReader()
        reader.SetFileName(inputfile)
        reader.Update()
        output1 = reader.GetOutput()
        # Specify the VTU output file name
        vtu_file = os.path.join(output_dir, filename)  #'red_airway-20361.vtu'
        # Create a VTU writer
        vtu_writer = vtk.vtkXMLUnstructuredGridWriter()
        vtu_writer.SetFileName(vtu_file)

        # Set the input as the output of the PVTU reader
        vtu_writer.SetInputData(output1)

        # Write the VTU file
        vtu_writer.Write()
    
    def find_peak(self, rm_unnecessary_pressure = True, time_over = -1 ,plot_ = False, plot_all_ = False):
        '''
        - Finding the peaks of PEEP and P-plateau
        - Input parameter:
            1. rm_unnecessary_pressure: ```bool```, it will determine whether the peaks which are close and in front of the begining of each pressure level remove or not.
            2. time_over : ```int```, only use when rm_unnecessary_pressure = True.  If the interval of two peaks are smaller than this value, then it the front peak will be ignored.
            3. plot_ : ```bool```, plot the pressure information only, without scaling.
            4. plot_all_: ```bool```, plot all the information, including flow and pressure, all the values are normalized. 
        - Output:
            1. ```np.array``` of P-plateau, where array[0] = index, array[1]= value
            2. ```np.array``` of PEEP, where array[0] = index, array[1]= value
        '''
        self.peak_pos = np.array([[]])
        self.peak_neg = np.array([[]])

        normal_result_P = self.pressure/np.max(self.pressure)
        normal_result_F = self.flow/np.max(self.flow)

        # Findind P- plateau
        peaks_neg, _ = find_peaks( -normal_result_F, distance=50, prominence=(0.5,None))
        peak_values_neg = normal_result_F[peaks_neg]
        peak_indices_neg =  self.time[peaks_neg]

        peak_start_index = np.array([])
        for peak_index in peaks_neg:
            difference = 0
            index = peak_index
            while difference <= 0 and index > 0:
                difference = normal_result_F[index] - normal_result_F[index -1]
                index = index -1
            if np.abs(normal_result_F[index]) >= 0.2:
                index = peak_index
                NotCross = True
                while NotCross and index > 0:
                    if normal_result_F[index] <= 0 and normal_result_F[index -1] >= 0:
                        NotCross = False
                    index = index -1
            peak_start_index = np.append(peak_start_index, index-1 )

        self.peak_pos = np.vstack((self.pressure_index[peak_start_index.astype(int)], self.pressure[peak_start_index.astype(int)]))


        # Findind PEEP
        peaks_pos, _ = find_peaks( normal_result_F, distance=50, prominence=(0.5,None))
        peak_values_pos = normal_result_F[peaks_pos]
        peak_indices_pos =  self.time[peaks_pos]

        peak_start_index0 = np.array([])
        for peak_index in peaks_pos:
            difference = 0
            index = peak_index
            while difference >= 0 and index > 0:
                difference = normal_result_F[index] - normal_result_F[index -1]
                index = index -1
            if np.abs(normal_result_F[index]) >= 0.2:
                index = peak_index
                NotCross = True
                while NotCross and index > 0:
                    if normal_result_F[index] >= 0 and normal_result_F[index -1] <= 0:
                        NotCross = False
                    index = index -1
            peak_start_index0 = np.append(peak_start_index0, index-1 )


        self.peak_neg = np.vstack((self.pressure_index[peak_start_index0.astype(int)], self.pressure[peak_start_index0.astype(int)]))

        # array for plotting
        time1 =  self.peak_pos[0,:]
        time2 =  self.peak_neg[0,:]
        all_time = np.append(time1,time2)
        all_time = np.sort(all_time)
        


        if rm_unnecessary_pressure:
            if time_over == -1:
                print("please define the time limit that if the time difference over it will be ignored")
            else:
                index_selected = np.array([round((all_time[0]-self.time[0])*self.frequency)])
                self.time_accepted = np.array([all_time[0]])
                for i in range(len(all_time)-1):
                    if all_time[i+1]-all_time[i] <= time_over:
                        self.time_accepted = np.append(self.time_accepted, all_time[i+1] )
                        index_selected = np.append(index_selected,round((all_time[i+1]-self.time[0])*self.frequency))
                
                pos_index = np.intersect1d(peak_start_index.astype(int), index_selected)
                neg_index = np.intersect1d(peak_start_index0.astype(int), index_selected)
                self.peak_pos = np.vstack((self.pressure_index[pos_index], self.pressure[pos_index]))
                self.peak_neg = np.vstack((self.pressure_index[neg_index], self.pressure[neg_index]))       
        else:
            self.time_accepted = all_time




        if plot_:
            plt.plot(self.time , self.pressure)
            if rm_unnecessary_pressure:
                plt.plot(self.pressure_index[pos_index], self.pressure[pos_index], 'ro')
                plt.plot(self.pressure_index[neg_index], self.pressure[neg_index], 'go')
            else:
                plt.plot(self.pressure_index[peak_start_index.astype(int)], self.pressure[peak_start_index.astype(int)], 'ro')
                plt.plot(self.pressure_index[peak_start_index0.astype(int)], self.pressure[peak_start_index0.astype(int)], 'go')   
                        
            plt.xlabel('Time')
            plt.ylabel('avg(pressure)')
            plt.title('Plot of CSV Column')
            plt.grid(True)
            plt.show()
        
        if plot_all_:
            plt.plot(self.time , normal_result_P,label='pressure (normalized)')
            plt.plot(self.time , normal_result_F,label='flow (normalized)')
            plt.plot(self.time[peak_start_index.astype(int)], normal_result_F[peak_start_index.astype(int)], 'ro',label='Plateau detected point')
            plt.plot(self.time[peak_start_index.astype(int)], normal_result_P[peak_start_index.astype(int)], 'bo',label='P plateau')
            plt.plot(self.time[peak_start_index0.astype(int)], normal_result_F[peak_start_index0.astype(int)], 'ro',label='PEEP detected point')
            plt.plot(self.time[peak_start_index0.astype(int)], normal_result_P[peak_start_index0.astype(int)], 'go',label='PEEP')
            if rm_unnecessary_pressure:
                plt.plot(self.time_accepted, normal_result_P[index_selected], 'y*',label='exported')
            plt.xlabel('Time')
            plt.ylabel('avg(pressure)')
            plt.title('Plot of CSV Column')
            plt.legend()
            plt.grid(True)
            plt.show()
        

        return self.peak_pos, self.peak_neg


               
    def indent(self,elem, level=0):
        '''
        - Funciton used for adjust the format of XML fils
        '''
        
        # Add indentation and line breaks to the XML tree
        indent_str = "\n" + level * "  "
        if len(elem):
            if not elem.text or not elem.text.strip():
                elem.text = indent_str + "  "
            if not elem.tail or not elem.tail.strip():
                elem.tail = indent_str
            for elem in elem:
                self.indent(elem, level + 1)
            if not elem.tail or not elem.tail.strip():
                elem.tail = indent_str
        else:
            if level and (not elem.tail or not elem.tail.strip()):
                elem.tail = indent_str

    def plot_output_config_range(self,dir,name,type,cut_point_1,cut_point_2):
        plt.plot(self.time , self.pressure)
        plt.plot(self.time_accepted , self.pressure[((self.time_accepted-self.time[0])*self.frequency).astype(int)],'ro',label = 'exported peaks')
        plt.axvline(x=cut_point_1, color='yellow', linestyle='--',label = 'start')
        plt.axvline(x=cut_point_2, color='green', linestyle='--',label = 'end')
        plt.xlabel('Time')
        plt.ylabel('avg(pressure)')
        plt.title('pressure level '+ str(name))
        plt.grid(True)
        plt.legend()
        file_name = "output_visualization_for_config_"+type+"_"+str(name)+".png"
        path = os.path.join(dir,file_name)
        plt.savefig(path)
        plt.close()
        


    def createXMLfile(self,inpur_dir, template_dir,template_name, output_dir,filename, shift = 0, createPVD = False, transfer2VTU = False, createConfig_sim = False, createConfig_meas = False,cut_sim = 150, cut_meas = 10, limit_sim = [0,-1],limit_meas = [0,-1]):
        '''
        - createXMLfile generate the PVD file included the group of pvtu files as the XML formate. 
        - it can also generate the vtu files for Holo-shed algorithm by making the input variable: transfer2VTU = True. 
        - Both the PVD file and VTU files will be stored in designated output directory.
        - input:
            1. ```str```: inpur_dir: directory that you store the original PVTU and VTU files
            2. ```str```: template_dir: directory of the config file's template
            3. ```str```: template_name: name of the config file's template
            4. ```str```: output_dir:  directory that you want to store the processed PVTU / VTU files
            5. ```str```: filename: PVD file_name , plese add the ".pvd" behind the name. e.g."output_SMART2.pvd"
            6. ```float```: shift: If you find the pressure extract from the PVD file (check from Paraview) are not located at your designated position, 
                             it might due to non-time-sync between your PVD files and your raw pressure data. Using this variable to shift the time (index) of each PVTU file,
                             and store the PVTU file into PVD file with correct time position. 
            7. ```bool```: transfer2VTU: In holo-shed alogorithm, the limitation is that it cannot read the PVTU file directly, therefore we have to first combine PVTU files of each 
                                   frame to one VTU file, then put the group of those pre-processed VTU file into "01_input folder" for holo-shed algorithm.
                                   Please check holo-shed library for the detail.
            8. ```bool```: createConfig: If true, the config files for holo-shed algorithm will be created, please put this .json file into "00_config" folder. notice that the json file is created by 
                     the sub-function: create_config_HoloShed_used(). If one want to modify the parameter regardless the timestep, one can modify them in the function directly.   
            9. ```list```: limit: use this define the upper and lower bound of the index, if limit = [0, -1], then it will loop over all the index.  Notice that for simulation case, the range is determined by frame index and for measurement case, it is determined by time (sec).
        - output:
            1. the ```.pvd``` file will be stored at designated output directory and naming as filename.
            2. if ```transfer2VTU = True```, then the ```.vtu``` files will also be stored at designated output directory, and naming as red_airway-*****.vtu, where ***** is the index value.
               Please notice that this option just for generate the ```.vtu``` files for holo-shed aldorithm, if you already have it, then you can de-activate it.
            3. if ```createConfig = True```, the config_file_*.json file will be stored at designated output directory. This is the neccessay for running the holo-shed algorithm.
        '''
        if not os.path.exists(os.path.join(output_dir)):
            os.makedirs(os.path.join(output_dir))

        if createPVD:
            datasets = []
            for time in self.time_accepted:
                if limit_sim[1] == -1:
                    if int(((time- shift)*self.frequency)-1) >= limit_sim[0] :
                        file = "test_SMART_E2_lung_new_V15000_awac_nlin-files/red_airway-"+ str(int(((time- shift)*self.frequency)-1))+".pvtu"
                        datasets.append({"timestep": str(time - shift), "group": "", "part": "0", "file": file})
                else:
                    if int(((time- shift)*self.frequency)-1) >= limit_sim[0] and int(((time- shift)*self.frequency)-1) <= limit_sim[-1]:
                        file = "test_SMART_E2_lung_new_V15000_awac_nlin-files/red_airway-"+ str(int(((time- shift)*self.frequency)-1))+".pvtu"
                        datasets.append({"timestep": str(time - shift), "group": "", "part": "0", "file": file})

            # Create the root element
            root = ET.Element("VTKFile")
            root.set("type", "Collection")
            root.set("version", "0.1")
            root.set("ByteOrder", "LittleEndian")

            # Create the Collection element
            collection = ET.SubElement(root, "Collection")

            for data in datasets:
                dataset = ET.SubElement(collection, "DataSet")
                for key, value in data.items():
                    dataset.set(key, value)

            # Create the XML tree
            tree = ET.ElementTree(root)

            # Set up indentation for pretty printing
            self.indent(root)
            
            # Write the XML tree to a file
            tree.write(os.path.join(output_dir, filename), encoding="utf-8", xml_declaration=True)

        if transfer2VTU:
            count = 0
            for time in self.time_accepted:

                if limit_sim[1] == -1:
                    if int(((time- shift)*self.frequency)-1) >=limit_sim[0] :
                        path_PVTU = os.path.join(inpur_dir,"red_airway-"+ str(int(((time- shift)*self.frequency)-1))+".pvtu")
                        self.generateVTUdata(path_PVTU ,output_dir ,"red_airway-"+ str(int(((time- shift)*self.frequency)-1))+".vtu")
                else:
                    if int(((time- shift)*self.frequency)-1) >= limit_sim[0] and int(((time- shift)*self.frequency)-1) <= limit_sim[-1]:
                        path_PVTU = os.path.join(inpur_dir,"red_airway-"+ str(int(((time- shift)*self.frequency)-1))+".pvtu")
                        self.generateVTUdata(path_PVTU ,output_dir ,"red_airway-"+ str(int(((time- shift)*self.frequency)-1))+".vtu")

                count = count+1
                if count % int(len(self.time_accepted)/20) == 0: 
                    print("already transfer",count, " / ",len(self.time_accepted))
        
        if createConfig_sim:
            file_index = 0
            time_accepted_index = (self.time_accepted-shift)*self.frequency-1
            if limit_sim[1] == -1:
                time_accepted_index = [x for x in time_accepted_index if 0 <= x <= time_accepted_index[-1]]
            else:
                time_accepted_index = [x for x in time_accepted_index if limit_sim[0]<= x <= limit_sim[1]]
            
            start_time = time_accepted_index[0]
            for i in range(len(time_accepted_index)):
                if i!=0:
                    interval = time_accepted_index[i] - time_accepted_index[i-1]
                    # the interval here is fixed at 150 frames currently, but one can modify if the interval of each level of Pressrue are different
                    if interval >= cut_sim:
                        end_time = time_accepted_index[i-1]
                        create_config_HoloShed_sim(int(start_time) ,int(end_time),template_dir ,template_name, output_dir ,"config_file_sim_"+str(file_index)+".json" )
                        self.plot_output_config_range(output_dir,file_index,'sim', start_time/self.frequency,end_time/self.frequency)
                        file_index = file_index +1
                        start_time = time_accepted_index[i]
            end_time = time_accepted_index[-1]
            create_config_HoloShed_sim(int(start_time) ,int(end_time),template_dir ,template_name, output_dir ,"config_file_sim_"+str(file_index)+".json" )
            self.plot_output_config_range(output_dir,file_index,'sim', start_time/self.frequency,end_time/self.frequency)
            print("generate Config files sucessfully")
        
        if createConfig_meas:
            index_list = np.array([])
            file_index = 0
            time_accepted_index = (self.time_accepted-shift)
            if limit_meas[1] == -1:
                time_accepted_index = [x for x in time_accepted_index if 0 <= x <= time_accepted_index[-1]]
            else:
                time_accepted_index = [x for x in time_accepted_index if limit_meas[0]<= x <= limit_meas[1]]
            
            for i in range(len(time_accepted_index)):
                if i!=0:
                    interval = time_accepted_index[i] - time_accepted_index[i-1]
                    # the interval here is fixed at 150 frames currently, but one can modify if the interval of each level of Pressrue are different
                    if interval <= cut_meas:
                        index_list = np.append(time_accepted_index[i-1], index_list)
                    else:
                        index_list = np.append(time_accepted_index[i-1], index_list)
                        index_list = index_list[::-1]
                        print("level_",str(file_index),": ",index_list)
                        create_config_HoloShed_meas(index_list, template_dir,template_name, output_dir ,"config_file_meas_"+str(file_index)+".json")

                        self.plot_output_config_range(output_dir,file_index,'meas', index_list[0],index_list[-1])
                        file_index = file_index +1
                        index_list = np.array([])
            index_list = np.append(time_accepted_index[i-1], index_list)
            index_list = index_list[::-1]
            print("level_",str(file_index),": ",index_list)
            create_config_HoloShed_meas(index_list, template_dir,template_name, output_dir ,"config_file_meas_"+str(file_index)+".json")

            self.plot_output_config_range(output_dir,file_index,'meas',index_list[0],index_list[-1])
            print("generate Config files sucessfully")

        


def calculate_decremental_pressure(peak ,time_thershold =10 , p_thershold = 100, NumValueCounted = -1, abandon_last = False, plot_ = False):
    '''
    - This funciton will create the list of average pressures along the whole pressure level. 
    - Input:

        1. ```np.array```, peak: input peak data, usually the upper or lower peak generated by ExtractPVDandPressure.find_peak(). 
        2. ```float``` , time_thershold: the time interval for each pressure level, can be used with p_thershold
        3. ```float``` , p_thershold: the pressure interval for each pressure level, can be used with time_thershold
        1. ```int```, NumValueCounted: numbers of data for computing the average. if ```NumValueCounted = -1```, then all the data will be used.
        2. ```bool```, abandon_last: if ```abandon_last - True```, then the last peak of each pressure level will be ignored.
        3. ```bool```, plot_: plot the pressure or not
    - Output:
        1. ```np.array``` of average pressure along whole pressure level 
    '''
    time = peak[0,:]
    pressure = peak[1,:]
    value = np.array([])
    value_array = np.array([])
    for i in range(len(pressure)-1):
        if abs(pressure[i+1] - pressure[i]) <= p_thershold  and  i != len(pressure)-2 and time[i+1]-time[i]<=time_thershold:
            value = np.append(value,pressure[i])
        else:
            if not abandon_last:
                value = np.append(value,pressure[i])
                
            if len(value) <= NumValueCounted:
                NumValueCounted = -1
                print("WARNING!,NumValueCounted is larger than the length of array, use all of the array for avergae instead")

            if NumValueCounted == -1:
                average = np.mean(value)
            else:
                array_counted = value[-NumValueCounted:]
                average = np.mean(array_counted)

            value_array = np.append( value_array, average)
            value = np.array([])
    if plot_:
        plt.plot(time , pressure)
        plt.xlabel('Time')
        plt.ylabel('avg(pressure)')
        plt.title('Plot of CSV Column')
        plt.grid(True)
        plt.show()
        

    return value_array
