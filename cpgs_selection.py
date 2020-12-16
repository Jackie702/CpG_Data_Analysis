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
        self.label_s = None
        
        self.discret_data = None
        self.info_gain_list = []
        self.filtered_ig_list = []
        self.qualified_feature_index = []
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
        self.label_s = pd.Series(self.label_data.iloc[:, 0])       # make label pandas.DataFrame info pandas.Series

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
        #self.print_data()
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

    # calculate a columns of info entropy + conditional entropy
    def calc_one_column_info_gain(self, data_column):
        p_x_list = data_column.value_counts() / data_column.shape[0]    # calculate the probility of each x element
        p_y_list = self.label_s.value_counts() / self.label_s.shape[0]          # calculate the probility of each y element

        print("#################################################################################################")
        print("\np_x_list:\n", p_x_list)
        print("\np_y_list:\n", p_y_list)

        H_x = 0
        for i in range(p_x_list.shape[0]):
            # print(f"p_x_list[{i}]={p_x_list[i]}")
            if p_x_list[i] == 0:
                continue
            else:
                H_x += -(p_x_list[i] * math.log(p_x_list[i], 2))

        print("\nH_x:\n", H_x)

        # 问题出在这里！！！！！！！！！！！！  data_column的index是TCGA-XX-XXXX-XX， 需要改为序号
        x_label_eql_0 = data_column[self.label_s[data_column.index] == 0]       # get numbers whose label is 0
        x_label_eql_1 = data_column[self.label_s[data_column.index] == 1]       # get numbers whose label is 1
        
        print(f"x_label_eql_0:\n {x_label_eql_0}")
        print(f"x_label_eql_1:\n {x_label_eql_1}\n")
        
        p_x_y_0 = x_label_eql_0.value_counts() / x_label_eql_0.shape[0]     # calculate p(X|Y=0)
        p_x_y_0 = np.array([n for n in p_x_y_0 if n != 0])                  # eliminate 0 elements for logarithm
        p_x_y_1 = x_label_eql_1.value_counts() / x_label_eql_1.shape[0]     # calculate p(X|Y=1)
        p_x_y_1 = np.array([n for n in p_x_y_1 if n != 0])                  # eliminate 0 elements for logarithm
        
        print(f"p_x_y_0:\n {p_x_y_0}")
        print(f"p_x_y_1:\n {p_x_y_1}\n")
        
        H_x_y_eql_0 = - (p_x_y_0 * np.log2(p_x_y_0)).sum()                  # calculate conditional entropy H(X|Y=0)
        H_x_y_eql_1 = - (p_x_y_1 * np.log2(p_x_y_1)).sum()                  # calculate conditional entropy H(X|Y=1)
        
        print("\nH_x_y_eql_0:", H_x_y_eql_0)
        print("\nH_x_y_eql_1:", H_x_y_eql_1)
        
        H_x_y = p_y_list[0] * H_x_y_eql_0 + p_y_list[1] * H_x_y_eql_1       # calculate conditon entropy H(X|Y=label)
        
        print("\nH_x_y:", H_x_y)
        print("#################################################################################################")

        return H_x + H_x_y
    
    # calculate infomation gain column by column, with sequence number in the list
    # [(index 0, ig 0), (index 1, ig 1), ..., (index N, ig N)]
    def calc_all_info_gain(self):
        for i in range(self.discret_data.shape[1]):
            ig = self.calc_one_column_info_gain(self.discret_data.iloc[:, i])
            self.info_gain_list.append((i, ig))
        return 

    def sort_and_filt(self):
        self.info_gain_list.sort(key=lambda x: x[1], reverse=True)  # sort the list in descending order
        
        # eliminate features whose info gain is less than threshold
        for i in range(len(self.filtered_ig_list)):
            if self.filtered_ig_list[i][1] >= self.threshold:
                self.filtered_ig_list.append(self.filtered_ig_list[i])
            else:
                break
        return

    # get combinations of feature by indices
    def get_feature_combs(self):
        self.qualified_feature_index = [x[0] for x in self.filtered_ig_list]    # get indices of qualified features
        index_combs = combinations(self.qualified_feature_index, self.subset_size)
        for i in range():
            
            for j in range(self.subset_size):
                    pass


if __name__ == "__main__":
    # create a argument parser
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument("-c", "--cpg", help="the file path of csv file with cpg data", default="./myNorm_0829model.csv")
    arg_parser.add_argument("-l", "--label", help="the file path of csv file with labels", default="./group.csv")
    arg_parser.add_argument("-t", "--threshold", help="the threshold for eliminating info gains", type=int, default=0.3)
    arg_parser.add_argument("-s", "--subset", help="the size of selected cpg subset", type=int, default=4)
    args = arg_parser.parse_args()
    
    # read and organize data
    cpg_data = data_container(args.cpg, args.label, args.threshold, args.subset)
    cpg_data.read_and_organize_data()
    cpg_data.save_transformed_csv()
    cpg_data.discretize_all()
    cpg_data.calc_all_info_gain()
    cpg_data.sort_and_filt()
    
    
    
    
    
    
    
    
    
    
    
        # ##### discretize raw data #####
    # def _discretize_data(self):
    #     self.discret_matrix = copy.deepcopy(self.cpg_data)
    #     intervals = [x/10 for x in range(11)]
    #     label_values = [y/10 for y in range(1, 11)]
    #     for i in range(self.discret_matrix.shape[1]):
    #         discretized_column = pd.cut(self.cpg_data.iloc[:, i], bins=intervals, labels=label_values)
    #         self.discret_matrix.iloc[:, i] = discretized_column
    #     return