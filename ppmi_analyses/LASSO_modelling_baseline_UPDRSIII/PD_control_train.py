# import needed packages
import numpy as np
import pandas as pd
from digipd_ml.supervised.classification import NestedCV

# set seed
np.random.seed(111)

# I/O
IN_DIR = "../data/"
OUTPUT_DIR = '../results/'
target = "UPDRS3_binary" # DIAGNOSIS
input_file = "data_cv_expr_4ML_" + target + ".tsv"
output_file = "results_nestedCV_" + target + ".csv"

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

sbatch pylauncher.sbatch PD_control_finalmodel_plot_PW.py CORUM mean linearSVM
sbatch pylauncher.sbatch PD_control_finalmodel_plot_PW.py CORUM median LogisticRegression
sbatch pylauncher.sbatch PD_control_finalmodel_plot_PW.py CORUM pca RandomForest
sbatch pylauncher.sbatch PD_control_finalmodel_plot_PW.py CORUM pathifier GradientBoosting
sbatch pylauncher.sbatch PD_control_finalmodel_plot_PW.py CORUM sd Adaboost

sbatch pylauncher.sbatch PD_control_finalmodel_plot_PW.py GOBP mean LogisticRegression
sbatch pylauncher.sbatch PD_control_finalmodel_plot_PW.py GOBP median RandomForest
sbatch pylauncher.sbatch PD_control_finalmodel_plot_PW.py GOBP pca RBFSVM
sbatch pylauncher.sbatch PD_control_finalmodel_plot_PW.py GOBP pathifier GradientBoosting
sbatch pylauncher.sbatch PD_control_finalmodel_plot_PW.py GOBP sd LogisticRegression

sbatch pylauncher.sbatch PD_control_finalmodel_plot_PW.py GOCC mean RBFSVM
sbatch pylauncher.sbatch PD_control_finalmodel_plot_PW.py GOCC median RBFSVM
sbatch pylauncher.sbatch PD_control_finalmodel_plot_PW.py GOCC pca LogisticRegression
sbatch pylauncher.sbatch PD_control_finalmodel_plot_PW.py GOCC pathifier linearSVM
sbatch pylauncher.sbatch PD_control_finalmodel_plot_PW.py GOCC sd GradientBoosting
