'''
Author: Ciyuan YU
Date: 2020-08-31 22:06:14
LastEditTime: 2020-09-23 22:58:33
LastEditors: Please set LastEditors
Description: Calculate CpG data to find the combination with largest impact on cervical tumor
FilePath: /CpG_data_analysis/cpgs_selection_v1.0.py
'''
import numpy as np
import pandas as pd
import math
from itertools import combinations
import copy
import argparse

class data_container:
    def __init__(self, cpg_file, label_file):
        self.cpg_path = cpg_file
        self.label_path = label_file
        self.cpg_obj = None
        self.label_obj = None
        
        self.cpg_sample_size = None
        self.cpg_features_size = None
        
        self.discret_matrix = None
        self.entropy_matrix = None

        self.feature_qualified = []
        
    
    ##### read in, organize and print data #####
    def _read_files(self):
        self.cpg_obj = pd.read_csv(self.cpg_path)
        self.label_obj = pd.read_csv(self.label_path)
        return
        
    def _organize_raw_data(self):
        # drop the sequential numbers of labels
        self.label_obj.drop(self.label_obj.columns[0], axis=1, inplace=True)
        self.label_obj.replace("Tumor", 1, inplace=True)
        self.label_obj.replace("Normal", 0, inplace=True)
        
        # convert the matrix of cpg data
        self.cpg_obj.set_index("Unnamed: 0", inplace=True)
        self.cpg_obj = pd.DataFrame(self.cpg_obj.values.T, 
                                    index=self.cpg_obj.columns, 
                                    columns=self.cpg_obj.index)
        
        # get shape info
        self.cpg_sample_size = self.cpg_obj.shape[0]
        self.cpg_features_size = self.cpg_obj.shape[1]
        return
    
    def print_data(self):
        print("\n####### cpg data #######\n")
        print(self.cpg_obj)
        print("\n####### labels #######\n")
        print(self.label_obj)
        print("\n####### Above are raw data we have after organization #######\n")
        return

    def read_and_organize_data(self):
        self._read_files()
        self._organize_raw_data()
        self.print_data()
        return
    
    def save_transformed_csv(self):
        self.cpg_obj.to_csv("cpg_data_transformed.csv")
        self.label_obj.to_csv("label_transformed.csv")
        return

    ##### discretize raw data #####
    def _discretize_data(self):
        self.discret_matrix = copy.deepcopy(self.cpg_obj)
        intervals = [x/10 for x in range(11)]
        label_values = [y/10 for y in range(1, 11)]
        for i in range(self.discret_matrix.shape[1]):
            discretized_column = pd.cut(self.cpg_obj.iloc[:, i], bins=intervals, labels=label_values)
            self.discret_matrix.iloc[:, i] = discretized_column
        return

    ##### compute entropy related information #####
    def _build_entropy_matrix(self):
        self.entropy_matrix = np.zeros((self.cpg_sample_size, self.cpg_features_size))
        return
    
    def _compute_one_sample_entropy(self, raw):
        # for i in range(self.cpg_features_size):
        #     self.entropy_matrix.iloc[raw, i] = 
        pass
    
    def _compute_all_entropy(self, raw):
        # for i in range(self.cpg_sample_size):
        pass   
            
        # return
    
    def _compute_features_permutations(self):
        pass
        


if __name__ == "__main__":
    # create a argument parser
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument("-c", "--cpg", help="the file path of csv file with cpg data", default="./myNorm_0829model.csv")
    arg_parser.add_argument("-l", "--label", help="the file path of csv file with labels", default="./group.csv")
    arg_parser.add_argument("-s", "--subset", help="the size of selected cpg subset", type=int, default=5)
    args = arg_parser.parse_args()
    
    # read and organize data
    cpg_data = data_container(args.cpg, args.label)
    cpg_data.read_and_organize_data()
    cpg_data.save_transformed_csv()
    
    