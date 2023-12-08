# import needed packages
import argparse
import numpy as np
import pandas as pd
from digipd_ml.supervised.classification import NestedCV

# set seed
np.random.seed(111)

# parse cmd arguments
# example python PD_control_train_PW.py CORUM mean
parser = argparse.ArgumentParser(prog = 'PD_control_train_PW',
                                 description = 'Nested CV at pathway level')
parser.add_argument('db', type=str, help='Aggregation database')
parser.add_argument('st', type=str, help='Aggregation statistic')
args = parser.parse_args()

# I/O
IN_DIR = "../data/"
OUTPUT_DIR = '../results/'
db = args.db
st = args.st
target = "UPDRS3_binary"
input_file = db + "_" + st + "_data_cv_expr_4ML_" + target + ".tsv"
output_file = db + "_" + st + "_results_nestedCV_" + target + ".csv"

# Reading the files
df = pd.read_table(IN_DIR + input_file, index_col=0)

# save diagnosis
y = df[target]
df = df.drop([target], axis=1)

# save features name and Patient_ID
features_name = [name for name in df.columns]
Patient_ID = list(df.index)

# Check that there is no more missing values
if np.sum(np.isnan(df)).sum() != 0:
    raise ValueError('There is missing data')

# transform data frame to numpy array
y = np.array(y)
X = np.array(df)

# perform nestedCV using class implemented in nestedcv.py file
# the balanced option will perform down sampling if set to False
# meaning that classes are not balanced
nestedCV = NestedCV()
nestedCV.fit(X, y, feat_select=True, normalize=True, balanced=True)

# save results into a csv
nestedCV.get_results().to_csv(OUTPUT_DIR + output_file)



