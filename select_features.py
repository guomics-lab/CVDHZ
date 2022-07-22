import pandas as pd
import numpy as np
import pdb
from scipy.stats import ttest_ind
from scipy.stats import zscore

file_feature = pd.read_csv('PBMC_protein.csv',index_col=0)
target_feature = pd.read_csv('PBMC_final_feature.csv',index_col=0)['protein']

temp_output = file_feature[target_feature]
thres = 0.5
low_NA_index = []
for i in temp_output.columns:
   temp_array = temp_output[i]
   #pdb.set_trace()
   if temp_array.isnull().sum() / len(temp_array) < thres:
       low_NA_index.append(i)

final_table = temp_output[low_NA_index]
final_table.fillna(final_table.min(axis=0), inplace=True)
final_table.to_csv('selected_PBMC_protein.csv')