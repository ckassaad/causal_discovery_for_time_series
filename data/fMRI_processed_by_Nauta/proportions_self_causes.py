from os import listdir
from os.path import isfile, join
import numpy as np

path_input = './returns/selected/'
files_input_name = [f for f in listdir(path_input) if isfile(join(path_input, f)) and not f.startswith('.')]
list_prop = []
for file_input_name in files_input_name:
    idx_ground_truth_file = file_input_name.split('timeseries')[1].split('.csv')[0]
    file_ground_truth_name = "sim" + idx_ground_truth_file + "_gt_processed"
    three_col_format_ground_truth = np.loadtxt(
        './ground_truth/' + file_ground_truth_name + '.csv',
        delimiter=',')
    count_self = 0
    for i in range(three_col_format_ground_truth.shape[0]):
        if int(three_col_format_ground_truth[i, 0]) == int(three_col_format_ground_truth[i, 1]):
            count_self = count_self + 1
    list_prop.append(str(count_self)+"/"+str(three_col_format_ground_truth.shape[0]))

print(list_prop)
# print(np.mean(list_prop))