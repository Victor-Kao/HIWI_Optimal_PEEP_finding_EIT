import csv
import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
from scipy.signal import find_peaks, peak_widths
import xml.etree.ElementTree as ET
from scipy import signal



class OptimalPEEP:
    def __init__(self,EIT_matrix, P_plateau, PEEP):
        self.EIT_matrix_ = EIT_matrix #np.array 2D
        self.P_plateau_ = P_plateau   #np.array 1D
        self.PEEP_ = PEEP             #np.array 1D
        self.PEEP_compliance_pixel_dict = {}#dict

        self.num_PEEP = self.EIT_matrix_.shape[0]
        self.num_row  = self.EIT_matrix_.shape[1]
        self.num_col  = self.EIT_matrix_.shape[2]

        self.cumulated_collpase = np.array([])
        self.cumulated_hyperdis = np.array([])

        self.percentage                 = np.zeros((self.num_PEEP,self.num_row ,self.num_col ))
        self.collapse_percentage        = np.zeros((self.num_PEEP,self.num_row ,self.num_col ))
        self.hyperdistention_percentage = np.zeros((self.num_PEEP,self.num_row ,self.num_col ))

        self.optimal_PEEP = -1
        self.optimal_value = -1

        #initialized indices
        self.indices = np.argwhere(self.collapse_percentage[0] == 1) # not meaningful, just some value never reached, so can obtain the empty array for initialized
        self.len_PEEP_flag = False
        self.len_P_plateau_flag = False
        self.initialize = False

        if self.num_PEEP == len(self.PEEP_):
            self.len_PEEP_flag = True
        else:
            print("WARNING! number of EIT_matrix != number of PEEP")
        
        if self.num_PEEP == len(self.P_plateau_):
            self.len_P_plateau_flag = True
        else:
            print("WARNING! number of EIT_matrix != number of P_plateau")

        if  self.len_PEEP_flag and self.len_P_plateau_flag:   
            print("initialized sucessfully")
            self.initialize = True
            

    def cal_compliance_pixel(self):
        self.compliance = np.zeros((self.num_PEEP, self.num_row , self.num_col ))
        for i in range(self.num_PEEP):
            for row in range(self.num_row):
                for col in range(self.num_col): 
                    if self.EIT_matrix_[i][row][col] ==0:
                        # give it some random number so it won't be 0, if not,sigularity might occur. 
                        self.EIT_matrix_[i][row][col] = 0.0000000000001*(np.random.random(1)+0.00001)
                    self.compliance[i][row][col] = self.EIT_matrix_[i][row][col] / (self.P_plateau_[i] - self.PEEP_[i])
        self.max_compliance_ = np.amax(self.compliance, axis=0)
        self.min_compliance_ = np.amin(self.compliance, axis=0)

        self.max_compliance = self.max_compliance_ - self.min_compliance_
        self.compliance = self.compliance - self.min_compliance_
        #print(self.max_compliance)

        print("calculate compliance sucessfully")
        pass


    def cal_percentage(self):
        for i in range(self.num_PEEP):
            self.percentage[i] = (self.max_compliance - self.compliance[i] ) *100 / self.max_compliance
            self.hyperdistention_percentage[i] = self.percentage[i]
            
            

            if i >= 1:
                if np.any(self.hyperdistention_percentage[i-1] == 0):
                    self.indices = np.concatenate((self.indices, np.argwhere(self.percentage[i-1] == 0)),axis=0)
                    
            if len(self.indices) != 0:
                #print(np.sort(self.indices[:, 1]))

                self.hyperdistention_percentage[i][self.indices[:, 0], self.indices[:, 1]] = 0
                self.collapse_percentage[i] = self.percentage[i] - self.hyperdistention_percentage[i]

        for i in range(self.num_PEEP):
            collapse = np.sum(self.collapse_percentage[i]*self.max_compliance ) / np.sum(self.max_compliance)
            #print(np.sum(self.collapse_percentage[i]*self.max_compliance ))
            hyperdis = np.sum(self.hyperdistention_percentage[i]*self.max_compliance ) / np.sum(self.max_compliance)
            #print(np.sum(self.hyperdistention_percentage[i]*self.max_compliance ))
            self.cumulated_collpase = np.append(self.cumulated_collpase , collapse)
            self.cumulated_hyperdis = np.append(self.cumulated_hyperdis , hyperdis)
        print("calculate percentage sucessfully")

    def find_optimal_PEEP(self):
        difference = self.cumulated_hyperdis - self.cumulated_collpase
        var_indices = np.where(np.logical_and(difference[:-1] >= 0, difference[1:] < 0))[0]
        if len(var_indices)==1:
            start_index = var_indices[0]
            end_index = var_indices[0]+1

            PEEP_start = self.PEEP_[start_index]
            PEEP_end = self.PEEP_[end_index]
            hyperdis_start = self.cumulated_hyperdis[start_index]
            hyperdis_end =  self.cumulated_hyperdis[end_index]
            collapse_start = self.cumulated_collpase[start_index]
            collapse_end = self.cumulated_collpase[end_index]

            m_h = ( hyperdis_end - hyperdis_start ) / (PEEP_end - PEEP_start)
            m_c = ( collapse_end - collapse_start ) / (PEEP_end - PEEP_start)

            self.optimal_PEEP = (collapse_start - hyperdis_start) / (m_h -  m_c) + PEEP_start
            self.optimal_value = m_h*(self.optimal_PEEP - PEEP_start) + hyperdis_start
            print("find optimal sucessfully")
            return self.optimal_PEEP
        else:
            print("WANRING! there is no intersection or more than one intersection ")
            print("please check cumulated_collpase and cumulated_hyperdis data ")
            print("find optimal failed")
            return 0


    def plot(self):
        # Create the plot
        plt.plot(self.PEEP_, self.cumulated_collpase  , label='Collpase')
        #print(self.cumulated_collpase)
        plt.plot(self.PEEP_, self.cumulated_hyperdis  , label='Hyperdistension')
        #print(self.cumulated_hyperdis)
        if self.optimal_PEEP >= 0 and self.optimal_value >= 0:
            plt.scatter(self.optimal_PEEP , self.optimal_value, label = 'optimal point', c= 'red')
            plt.text(self.optimal_PEEP, self.optimal_value, f'({self.optimal_PEEP:.2f}, {self.optimal_value:.2f})')
        # Add labels and title
        plt.xlabel('PEEP (cmH20)')
        plt.ylabel('Cumulated percentage (%)')
        plt.title('Estimation plot by EIT')
        # Display the plot
        plt.legend()
        plt.show()
     
    def analysis(self, plot = False):
        self.cal_compliance_pixel()
        self.cal_percentage()
        self.find_optimal_PEEP()
        if plot:
            self.plot()



def EITimport(input_file):
    '''
    Use funcntion to import the EIT matrix. 
    '''
    loaded_data = np.loadtxt(input_file,delimiter=',')
    n_arrays_import, flattened_shape = loaded_data.shape
    import_EIT_ = loaded_data.reshape((n_arrays_import,flattened_shape,1))
    import_EIT_ = import_EIT_[:,~np.isnan(import_EIT_).any(axis=0)] 
    shape_0_,shape_1_ = import_EIT_.shape
    return import_EIT_,shape_0_,shape_1_


def findOptimal(file_dict_, P_pleateau_, PEEP_, peak_type_ = 'single', selected_peak_ = -1, plot_ = True):
    '''
    - Input:
        1. ```dict```, file_dict_: please check the SMART2_case to see how to construct the dict with correct format
        2. ```np.array```, P_pleateau_: array of P pleateau along whole pressure level. It can be genertated by calculate_decremental_pressure() from PEEP_Ppleateau_detect_lib
        3. ```np.array```, PEEP_: array of PEEP along whole pressure level. It can be genertated by calculate_decremental_pressure() from PEEP_Ppleateau_detect_lib
        4. ```str```, frame_type_ : If frame type is ```'single'```, then the algorithm will only choose one peak from each pressure level for further computing. 
                                   If frame type is ```'avg'```, then it will average the numbers of EIT data for further computing  
        5. ```int```, selected_frame_: If frame_type_ = ```'single'```, then the selected_frame_ will determine which Peak in each pressure level will be selected for computing.
                                       Notice that it is counted from the end, so selected_frame_ = -1 means that it will select the last peak for computing. 
                                       If frame_type_ = ```'avg'```, this parameter will be ignored.
        6. ```bool```, plot_: plot the result 
    '''

    #check the size
    dict_length = len(file_dict_)
    P_pleateau_length = len(P_pleateau_)
    PEEP_length = len(PEEP_)
    #initialize and 
    shape_0_size_array = np.array([])
    shape_1_size_array = np.array([])
    EIT_matrix = np.array([])

    #import the data
    if dict_length == P_pleateau_length and P_pleateau_length == PEEP_length:

        
        for key in file_dict_:
            path_ = file_dict_[key]
            _ , shape_0, shape_1 = EITimport(path_)
            shape_0_size_array = np.append(shape_0_size_array,shape_0)
            shape_1_size_array = np.append(shape_1_size_array,shape_1)

        EIT_P_pleateau_ = np.zeros([len(shape_1_size_array),1,shape_1])
        EIT_PEEP_       = np.zeros([len(shape_1_size_array),1,shape_1])

        if peak_type_ == 'single':
            if all(num == shape_1_size_array[0] for num in shape_1_size_array):
                Num_peak = np.min(shape_0_size_array)
                #initalize the matrix for storing EIT value
                #flatten the 2D matrix to 1D array previosuly 
                #EIT_P_pleateau_ = np.zeros([len(shape_1_size_array),1,shape_1])
                #EIT_PEEP_       = np.zeros([len(shape_1_size_array),1,shape_1])
                #import the value for selected frame
                if abs(selected_peak_) >= round(Num_peak/2)-1:
                    selected_peak_ =  round(Num_peak/2)-1
                    print("WARNING! selected frame is over the range, select the maximum acceptable frame = ",round(Num_peak/2)-1," instead")
                frame_pos_ = 2*selected_peak_+1
                frame_neg_ = frame_pos_-1

                num_counted = 0
                for key in file_dict_:
                    path_ = file_dict_[key]
                    EIT_value ,_ , _ = EITimport(path_)
                    EIT_P_pleateau_[num_counted] = EIT_value[frame_pos_][:][:]
                    EIT_PEEP_[num_counted] = EIT_value[frame_neg_][:][:]
                    num_counted = num_counted +1
                print("import the EIT value sucessfully")
                EITdiff_ = EIT_P_pleateau_ - EIT_PEEP_
            else:
                print("ERROR! the size of EIT matrix are not match ")

        elif peak_type_ == 'avg':
                num_counted = 0
                for key in file_dict_:
                    path_ = file_dict_[key]
                    EIT_value ,_ , _ = EITimport(path_)
                    #print(EIT_value.shape)
                    Peak_neg = EIT_value[0::2][:][:]
                    Peak_pos = EIT_value[1::2][:][:]
                    EIT_P_pleateau_[num_counted] = np.mean(Peak_pos,axis =0 )
                    EIT_PEEP_[num_counted] = np.mean(Peak_neg,axis =0 )
                    num_counted = num_counted +1
                print("import the EIT value sucessfully")

                #for i in range(len(EIT_P_pleateau_)):
                #    minimum_EIT_P_pleateau_ = np.min(EIT_P_pleateau_)
                #    minimum_EIT_PEEP_ = np.min(EIT_PEEP_)
                #    if minimum_EIT_P_pleateau_ <= 0:
                #        EIT_P_pleateau_[i] = EIT_P_pleateau_[i]- minimum_EIT_P_pleateau_
                #    if minimum_EIT_PEEP_ <= 0:
                #        EIT_PEEP_[i] = EIT_PEEP_[i]- minimum_EIT_PEEP_
                    
                
                EITdiff_ = EIT_P_pleateau_ - EIT_PEEP_
                #print(EIT_P_pleateau_)


        else: 
            print("ERROR!! Please input the correct peak_type: 'single' or 'avg'.")
        
        print("P_pleateau = ", P_pleateau_/100)
        print("PEEP = ", PEEP_/100)
        patient1 = OptimalPEEP(EITdiff_,P_pleateau_/100, PEEP_/100)
        patient1.analysis(plot= plot_)

    else:
        print("ERROR! the size of the files are not matched the size of Pressure list")
        pass