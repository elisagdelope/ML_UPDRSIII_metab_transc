# import needed packages
import numpy as np
import pandas as pd
import shap
from digipd_ml.supervised.classification import NestedCV
from digipd_ml.plot_utils import plot_shap_values
from digipd_ml.plot_utils import plot_from_shap_values
from digipd_ml.utils import feature_importances_shap_values
from digipd_ml.utils import feature_importances_from_shap_values
from sklearn import metrics
import argparse

# set seed
np.random.seed(111)

# parse cmd arguments
# example python PD_control_finalmodel_plot_notest.py lm-time LogisticRegression
parser = argparse.ArgumentParser(prog='PD_control_finalmodel_plot_notest.py',
                                 description='Nested CV at pathway level')
parser.add_argument('long_feat', type=str, help='Longitudinal feature')
parser.add_argument('model', type=str, help='Name of final model')
args = parser.parse_args()

# I/O
INPUT_DIR = "../data/"
OUTPUT_DIR = '../results/'
target = "UPDRS3_binary"  # DIAGNOSIS
model_name = args.model
temp_ft = args.long_feat
input_data_file = "data_" + temp_ft + "_expr_4ML_" + target + ".tsv"
out_ft_file = temp_ft + "_top_features_shap_" + target + ".csv"
out_shap_plot = temp_ft + "_shap_plot_" + target + ".pdf"

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

# shap values analysis
if model_name not in ["RBFSVM", "Adaboost"]:
    print("Shaps from explainer")

    # shap plot
    plot_shap_values(model['model'], X_processed, save_fig=True,
                     plot_title="Gene level features",
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
                          plot_title="Gene level features",
                          name_file=out_shap_plot, path=OUTPUT_DIR)

    # feature importance
    top_features = feature_importances_from_shap_values(shap_values, X_processed, n=20)
top_features.to_csv(OUTPUT_DIR + out_ft_file, index=False)

print("end")



sbatch pylauncher.sbatch PD_control_finalmodel_plot_notest.py sd RBFSVM
sbatch pylauncher.sbatch PD_control_finalmodel_plot_notest.py lm-time LogisticRegression

sbatch pylauncher.sbatch PD_control_finalmodel_plot_notest_PW.py GOBP mean sd linearSVM
sbatch pylauncher.sbatch PD_control_finalmodel_plot_notest_PW.py GOBP mean lm-time LogisticRegression
sbatch pylauncher.sbatch PD_control_finalmodel_plot_notest_PW.py GOBP median sd Adaboost
sbatch pylauncher.sbatch PD_control_finalmodel_plot_notest_PW.py GOBP median lm-time Adaboost
sbatch pylauncher.sbatch PD_control_finalmodel_plot_notest_PW.py GOBP pca sd linearSVM
sbatch pylauncher.sbatch PD_control_finalmodel_plot_notest_PW.py GOBP pca lm-time LogisticRegression
sbatch pylauncher.sbatch PD_control_finalmodel_plot_notest_PW.py GOBP pathifier sd LogisticRegression
sbatch pylauncher.sbatch PD_control_finalmodel_plot_notest_PW.py GOBP pathifier lm-time RBFSVM
sbatch pylauncher.sbatch PD_control_finalmodel_plot_notest_PW.py GOBP sd sd RBFSVM
sbatch pylauncher.sbatch PD_control_finalmodel_plot_notest_PW.py GOBP sd lm-time LogisticRegression

sbatch pylauncher.sbatch PD_control_finalmodel_plot_notest_PW.py GOCC mean sd LogisticRegression
sbatch pylauncher.sbatch PD_control_finalmodel_plot_notest_PW.py GOCC mean lm-time LogisticRegression
sbatch pylauncher.sbatch PD_control_finalmodel_plot_notest_PW.py GOCC median sd LogisticRegression
sbatch pylauncher.sbatch PD_control_finalmodel_plot_notest_PW.py GOCC median lm-time linearSVM
sbatch pylauncher.sbatch PD_control_finalmodel_plot_notest_PW.py GOCC pca sd RBFSVM
sbatch pylauncher.sbatch PD_control_finalmodel_plot_notest_PW.py GOCC pca lm-time Adaboost
sbatch pylauncher.sbatch PD_control_finalmodel_plot_notest_PW.py GOCC pathifier sd Adaboost
sbatch pylauncher.sbatch PD_control_finalmodel_plot_notest_PW.py GOCC pathifier lm-time RandomForest
sbatch pylauncher.sbatch PD_control_finalmodel_plot_notest_PW.py GOCC sd sd GradientBoosting
sbatch pylauncher.sbatch PD_control_finalmodel_plot_notest_PW.py GOCC sd lm-time LogisticRegression

sbatch pylauncher.sbatch PD_control_finalmodel_plot_notest_PW.py CORUM mean sd RBFSVM
sbatch pylauncher.sbatch PD_control_finalmodel_plot_notest_PW.py CORUM mean lm-time LogisticRegression
sbatch pylauncher.sbatch PD_control_finalmodel_plot_notest_PW.py CORUM median sd Adaboost
sbatch pylauncher.sbatch PD_control_finalmodel_plot_notest_PW.py CORUM median lm-time LinearSVM
sbatch pylauncher.sbatch PD_control_finalmodel_plot_notest_PW.py CORUM pca sd GradientBoosting
sbatch pylauncher.sbatch PD_control_finalmodel_plot_notest_PW.py CORUM pca lm-time Adaboost
sbatch pylauncher.sbatch PD_control_finalmodel_plot_notest_PW.py CORUM pathifier sd RBFSVM
sbatch pylauncher.sbatch PD_control_finalmodel_plot_notest_PW.py CORUM pathifier lm-time RBFSVM
sbatch pylauncher.sbatch PD_control_finalmodel_plot_notest_PW.py CORUM sd sd LogisticRegression
sbatch pylauncher.sbatch PD_control_finalmodel_plot_notest_PW.py CORUM sd lm-time Adaboost
