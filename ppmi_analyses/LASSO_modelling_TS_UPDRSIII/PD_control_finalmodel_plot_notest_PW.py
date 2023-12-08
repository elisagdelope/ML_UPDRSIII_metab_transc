# import needed packages
import argparse
import numpy as np
import pandas as pd
import shap
from digipd_ml.supervised.classification import NestedCV
from digipd_ml.plot_utils import plot_shap_values
from digipd_ml.plot_utils import plot_from_shap_values
from digipd_ml.utils import feature_importances_shap_values
from digipd_ml.utils import feature_importances_from_shap_values
from sklearn import metrics

# set seed
np.random.seed(111)

# parse cmd arguments
# example python PD_control_finalmodel_plot_notest_PW.py CORUM mean lm-time LogisticRegression
parser = argparse.ArgumentParser(prog='PD_control_finalmodel_plot_notest_PW.py',
                                 description='Nested CV at pathway level')
parser.add_argument('db', type=str, help='Aggregation database')
parser.add_argument('st', type=str, help='Aggregation statistic')
parser.add_argument('long_feat', type=str, help='Longitudinal feature')
parser.add_argument('model', type=str, help='Name of final model')
args = parser.parse_args()

# I/O
INPUT_DIR = "../data/"
OUTPUT_DIR = '../results/'
model_name = args.model
db = args.db
st = args.st
temp_ft = args.long_feat
target = "UPDRS3_binary"   # DIAGNOSIS
input_data_file = db + "_" + st + "_" + temp_ft + "_data_expr_4ML_" + target + ".tsv"
out_ft_file = db + "_" + st + "_" + temp_ft + "_top_features_shap_" + target + ".csv"
out_shap_plot = db + "_" + st + "_" + temp_ft + "_shap_plot_" + target + ".pdf"


# Reading the files
df = pd.read_table(INPUT_DIR + input_data_file, index_col=0)

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

# fit single best model
nestedCV = NestedCV()
model = nestedCV.fit_single_model(
    X, y, name_model=model_name, normalize=True, feat_select=True, balanced=True)

# shap values analysis
X_processed = model['selecter'].transform(model['scaler'].transform(X))  # apply scaler and selecter
feat_names = model['selecter'].get_feature_names_out(features_name)  # features selected (index ordered)
X_processed = pd.DataFrame(data=X_processed, columns=feat_names)

if "GOBP" in X_processed.columns[0]:
    pw_names = [string.replace("GOBP_", "") for string in X_processed.columns]
    pw_names = [string.replace("_", " ") for string in pw_names]
    pw_names = [string.capitalize() for string in pw_names]
elif "GOCC" in X_processed.columns[0]:
    pw_names = [string.replace("GOCC_", "") for string in X_processed.columns]
    pw_names = [string.replace("_", " ") for string in pw_names]
    pw_names = [string.capitalize() for string in pw_names]
else:
    pw_names = None

# shap values analysis
if model_name not in ["RBFSVM", "Adaboost"]:
    print("Shaps from explainer")

    # shap plot
    plot_shap_values(model['model'], X_processed, save_fig=True, names_list=pw_names,
                     plot_title="Pathway level " + st + " features",
                     name_file=out_shap_plot, path=OUTPUT_DIR)

    # feature importance
    top_features = feature_importances_shap_values(model['model'], X_processed, n=20)

else:
    print("Shaps from KernelExplainer")
    # generate shap values
    explainer = shap.KernelExplainer(model['model'].predict, shap.kmeans(X_processed, 100))
    shap_values = explainer.shap_values(X_processed)

    # shap plot
    plot_from_shap_values(shap_values, X_processed, save_fig=True,
                          plot_title="Pathway level " + st + " features", names_list=pw_names,
                          name_file=out_shap_plot, path=OUTPUT_DIR)

    # feature importance
    top_features = feature_importances_from_shap_values(shap_values, X_processed, n=20)
top_features.to_csv(OUTPUT_DIR + out_ft_file, index=False)

print("end")
