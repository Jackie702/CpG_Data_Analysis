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

# the switch of printing logs
test_open = 1

class data_container:
    def __init__(self, cpg_file = "./myNorm_0829model.csv", 
                 label_file = "./group.csv", 
                 threshold = 0.3, 
                 subset_num = 4):
        self.cpg_path = cpg_file
        self.label_path = label_file
        self.ig_threshold = threshold
        self.subset_size = subset_num
        
        self.cpg_data = None
        self.label_data = None
        self.np_data = None
        self.np_label = None
        
        self.discret_data = None
        self.label_entropy = None
        self.info_gain_list = []
        self.filtered_ig_list = []
        self.qualified_feature_index = []
        self.candidate_feature_dict = {}
        self.feature_combs = []

        self.feature_qualified = []
        
    
    ##### read in, organize and print data #####
    def _read_files(self):
        self.cpg_data = pd.read_csv(self.cpg_path)
        self.label_data = pd.read_csv(self.label_path)
        return
        
    def _organize_raw_data(self):
        # drop the sequential numbers of labels
        self.label_data.drop(self.label_data.columns[0], axis=1, inplace=True)
        self.label_data.replace("Tumor", 1, inplace=True)
        self.label_data.replace("Normal", 0, inplace=True)
        self.label_data = pd.Series(self.label_data.iloc[:, 0])       # make label pandas.DataFrame info pandas.Series

        # convert the matrix of cpg data
        self.cpg_data.set_index("Unnamed: 0", inplace=True)
        self.cpg_data = pd.DataFrame(self.cpg_data.values.T, 
                                    index=self.cpg_data.columns, 
                                    columns=self.cpg_data.index)
        self.cpg_data = self.cpg_data.reset_index(drop = True)      # reset cpg data index from TCGA-XX-XXXX-XX to 0 ~ n
        return
    
    def print_data(self):
        print("\n####### cpg data #######\n")
        print(self.cpg_data)
        print("\n####### labels #######\n")
        print(self.label_data)
        print("\n####### Above are raw data we have after organization #######\n")
        return

    def read_and_organize_data(self):
        self._read_files()
        self._organize_raw_data()
        
        if test_open:
            self.print_data()
        return
    
    def save_transformed_csv(self):
        self.cpg_data.to_csv("cpg_data_transformed.csv")
        self.label_data.to_csv("label_transformed.csv")
        return

    # discretize all raw data
    def discretize_all(self):
        self.discret_data = copy.deepcopy(self.cpg_data)
        intervals = [x/10 for x in range(11)]
        disc_values = [y for y in range(1, 11)]
        for i in range(self.cpg_data.shape[1]):
            disc_column = pd.cut(self.cpg_data.iloc[:, i], bins=intervals, labels=disc_values)
            self.discret_data.iloc[:, i] = disc_column
        return
    
    def calc_entropy(self, data):
        prob = pd.value_counts(data) / data.shape[0]
        return - sum(prob * np.log2(prob))
                
    def calc_info_gain(self, data, feature, condition):
        e_feature = data.groupby(feature).apply(lambda x: self.calc_entropy(x[condition]))
        p_data = pd.value_counts(data[feature]) / data[feature].shape[0]
        e_condition = sum(e_feature * p_data)
        print(f"e_feature:\n{e_feature}\np_data:\n{p_data}\ne_condition:\{e_condition}\n")
        return self.calc_entropy(data[condition]) - e_condition
    
    def calc_all_feature_info_gain(self):
        self.discret_data.insert(self.discret_data.shape[1], "label", self.label_data)
        print(self.discret_data)
        feature_set = set(self.discret_data.columns[:-1])
        for feature in feature_set:
            print(f"\n\n######## {feature} start ##########")
            info_gain = self.calc_info_gain(self.discret_data, feature, "label")
            print(f"feature_name:\n{feature}\ninfo_gain:\n{info_gain}\n")
            self.info_gain_list.append(info_gain)
        return
    
    def filt_and_sort(self):
        for i in range(len(self.info_gain_list)):
            ig = self.info_gain_list[i]
            if ig > self.ig_threshold:
                self.filtered_ig_list.append((i, ig))
        self.filtered_ig_list.sort(key=lambda x: x[1], reverse=True)
    
    def ig_list_to_dict(self):
        for item in self.filtered_ig_list:
            self.candidate_feature_dict[item[0]] = item[1]
    
    def comb_feature_by_SFFS_algo(self):
        for item in self.filtered_ig_list:
            pass
    
    def cpg_test(self):
        self.read_and_organize_data()
        self.discretize_all()
        self.calc_all_feature_info_gain()
        self.filt_and_sort()
        self.ig_list_to_dict()


if __name__ == "__main__":
    # create a argument parser
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument("-c", "--cpg", help="the file path of csv file with cpg data", default="./myNorm_0829model.csv")
    arg_parser.add_argument("-l", "--label", help="the file path of csv file with labels", default="./group.csv")
    arg_parser.add_argument("-t", "--threshold", help="the threshold for eliminating info gains", type=int, default=0.3)
    arg_parser.add_argument("-s", "--subset", help="the size of selected cpg subset", type=int, default=4)
    args = arg_parser.parse_args()
    
    cpg = data_container()
    cpg.cpg_test()





    # def calc_label_entropy(self):
    #     self.label_entropy = self.label_data.value_counts() / self.label_data.shape[0]
                
    # def calc_one_column_info_gain(self, index):
    #     np_H_y_x = np.zeros((1, 10))
    #     np_comp_mat = np.insert(self.np_data[:, index], 1, values=self.np_label, axis=1)
        
    #     element_eql_0 = [0] * 10
    #     element_eql_1 = [0] * 10
    #     p_y_x = [0] * 10
    #     for i in range(1, 11):
    #         for j in range(np_comp_mat.shape[0]):
    #             if np_comp_mat[j, 0] == i and np_comp_mat[j, 1] == 0:
    #                 element_eql_0[i - 1] += 1
    #             else:
    #                 element_eql_1[i - 1] += 1
        
    #     for k in range(10):
    #         p_y_x[k] = -(element_eql_0[k] * math.log(element_eql_0[k], 2) + 
    #                       element_eql_1[k] * math.log(element_eql_1[k], 2))
        
    #     for l in range(10):
            
        
    
    # def calc_all_info_gain(self):
    #     self.np_label = self.label_data.values        # get np array form of label
    #     self.np_data = self.discret_data.values       # get np array form of data
        
    #     for i in range(self.np_data.shape[1]):
    #         self.info_gain_list.append((i, self.calc_one_column_info_gain(i)))