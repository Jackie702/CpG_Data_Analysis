'''
Author: Ciyuan Yu
Date: 2020-09-23 21:31:12
LastEditTime: 2020-10-15 10:07:15
LastEditors: Please set LastEditors
Description: Use to test functinons in cpgs_selection
FilePath: /CpG_data_analysis/function_test.py
'''
import numpy as np
import pandas as pd
import math
from itertools import combinations
import copy
import argparse

titles = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i']
raw_data = [
[0.4221829,	    0.321612292,	0.665429392,	0.652345153,	0.25949937,	    0.256704268,	0.386967264,	0.13483489,	    0.417532239],
[0.945748086,	0.894798944,    0.938502461,	0.851027982,	0.594604807,	0.815685268,	0.825218378,	0.810564628,	0.879533024],
[0.864087523,	0.811501156,	0.970861069,	0.903106497,	0.504731955,	0.685082918,	0.713167475,	0.743161052,	0.779790625],
[0.061871919,	0.46760436,	    0.713360144,	0.669695774,	0.469950971,	0.05044445,	    0.067819408,	0.130357035,	0.373938953],
[0.796113932,	0.683169729,	0.908077088,	0.726079901,	0.547898799,	0.541301129,	0.716850749,	0.768435484,	0.745048613],
[0.874061301,	0.875764423,	0.976813773,	0.876217673,	0.566057263,	0.803680416,	0.779781403,	0.799567009,	0.798618217],
[0.707149151,	0.490108873,	0.780333604,	0.57769059,	    0.375718577,	0.284511938,	0.58896177,	    0.57597191,	    0.557430542],
[0.744540349,	0.292744219,	0.853421036,	0.588890036,	0.482925415,	0.473101265,	0.663032284,	0.201004742,	0.433609226],
[0.600824257,	0.406593235,	0.804926209,	0.53487398,	    0.223035994,	0.209517505,	0.427012183,	0.046816097,	0.38845731],
[0.792151961,	0.62859938,	    0.934155413,	0.751843078,	0.567949285,	0.616704237,	0.686977885,	0.714003439,	0.765023829],
[0.884533373,	0.80845848,	    0.907055287,	0.843609766,    0.588229675,	0.712064277,	0.775988471,	0.799285689,	0.732402697],
[0.838003757,	0.671983705,	0.852500265,	0.677227726,	0.520235863,	0.536002228,	0.687905423,	0.67037494,	    0.676818259],
[0.955328929,	0.913629428,	0.952332522,	0.859870423,	0.616041776,	0.80491712,	    0.789882716,	0.960767907,	0.889936348],
[0.592801,	    0.476418061,	0.969831435,	0.787532451,	0.549660282,	0.360043341,	0.659717221,	0.698255048,	0.503217119],
[0.86888519,	0.794145655,	0.964330332,	0.808215488,	0.544354174,	0.732794098,	0.762759475,	0.718958576,	0.795738141]
]
data_df = pd.DataFrame(raw_data)
data_df.columns = titles

labels = [1, 0, 1, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 1]
labels_s = pd.Series(labels)

def discret_column(df):
    df_new = copy.deepcopy(df)
    intervals = [x/10 for x in range(11)]
    disc_values = [y for y in range(1, 11)]
    for i in range(df.shape[1]):
         disc_column = pd.cut(df.iloc[:, i], bins=intervals, labels=disc_values)
         df_new.iloc[:, i] = disc_column
    return df_new

def calc_one_column_info_gain(data_column, test_print=True):
    ig = copy.deepcopy(data_column)
    p_x_list = ig.value_counts() / ig.shape[0]
    p_y_list = labels_s.value_counts() / labels_s.shape[0]
    p_x_y_0 = data_column[labels_s[:]==0].value_counts() / labels_s.shape[0]
    p_x_y_1 = data_column[labels_s[:]==1].value_counts() / labels_s.shape[0]

    if test_print:
    #     print("\n***** Probs of data cells：*****")
    #     print("p_x\tp_y\tp_x_y_0\tp_x_y_1")
    #     for i in range(ig.shape[0]):
    #         if(len(p_x_list)-1 >= i):  print("%.2f" % p_x_list[i], end="\t")
    #         if(len(p_y_list)-1 >= i):  print("%.2f" % p_y_list[i], end="\t")
    #         if(len(p_x_y_0)-1 >= i):   print("%.2f" % p_x_y_0[i], end="\t")
    #         if(len(p_x_y_1)-1 >= i):   print("%.2f" % p_x_y_1[i])

        print("\np_x_list:\n", p_x_list)
        print("\np_y_list:\n", p_y_list)
        print("\np_x_y_eql_0:\n", p_x_y_0)
        print("\np_x_y_eql_1:\n", p_x_y_1)

    for i in range(ig.shape[0]):
        entropy_x = 0
        entropy_x_y = 0

        # calculate the single entropy h(x)
        x = ig[i]
        p_x = p_x_list[x]
        entropy_x = -p_x * math.log(p_x, 2)

        # calculate the h(x|y)
        y = labels_s[i]
        p_y = p_y_list[y]
        factor = p_y
        if y == 0:
            factor *= p_x_y_0[x]
        elif y == 1:
            factor *= p_x_y_1[x]
        else:
            raise ValueError(f"Wrong label value in line {i+1}")
        entropy_x_y = -factor * math.log(factor, 2)

        ig[i] = entropy_x + entropy_x_y

    print("\n******** Ig of one column *********")
    print(ig)

    #return ig
        

# def calc_all_info_gain(disc_df, labels_s):
#     info_gain_mat = copy.deepcopy(disc_df)
#     for i in range(disc_df.shape[1]):
#         value_counts = disc_df.iloc[:, i].value_counts()
#         for j in range(disc_df.shape[0]):
#             p_x = disc_df.iloc[i, j]
#             info_gain_mat.iloc[i, j] = 
#     pass

if __name__ == "__main__":
    df_new = discret_column(data_df)
    print(f"\n********* Discretized matrix ***********\n{df_new}\n")
    calc_one_column_info_gain(df_new.iloc[:, 1])
    