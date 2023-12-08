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
INPUT_DIR = "../../../data/ts/classification/"
OUTPUT_DIR = '../../../results/ts/classification/'
model_name = args.model
temp_ft = args.long_feat
input_data_file = "data_" + temp_ft + "_metab_4ML_UPDRS__3_binary.tsv"
out_ft_file = temp_ft + "_top_features_shap_UPDRS__3_binary.csv"
out_shap_plot = temp_ft + "_shap_plot_UPDRS__3_binary.pdf"

# Reading the files
df = pd.read_table(INPUT_DIR + input_data_file, index_col=0)

# save diagnosis
y = df['UPDRS__3_binary']
df = df.drop(['UPDRS__3_binary'], axis=1)

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
    X, y, name_model=model_name, normalize=True, feat_select=True, balanced=False)

# annotation for metabolite-level data
annotation_file = "../../../../data/00-cleansing/chemical_annotation.tsv"
annotation_df = pd.read_table(annotation_file)
annotation_df = annotation_df.set_index("ANALYSIS_ID")

# shap values analysis
X_processed = model['selecter'].transform(model['scaler'].transform(X))  # apply scaler and selecter
feat_names = model['selecter'].get_feature_names_out(features_name)  # features selected (index ordered)
X_processed = pd.DataFrame(data=X_processed, columns=feat_names)
chem_names_list = list(annotation_df.loc[feat_names]["CHEMICAL_NAME"])

if model_name not in ["RBFSVM", "Adaboost"]:
    print("Shaps from explainer")

    # shap plot
    plot_shap_values(model['model'], X_processed, save_fig=True, names_list=chem_names_list,
                     plot_title="Metabolite level features",
                     name_file=out_shap_plot, path=OUTPUT_DIR)

    # feature importance
    top_features = feature_importances_shap_values(model['model'], X_processed, names_list=chem_names_list, n=20)

else:
    print("Shaps from KernelExplainer")
    # generate shap values
    if X_processed.shape[0] > 100:
        explainer = shap.KernelExplainer(model['model'].predict, shap.kmeans(X_processed, 100))
    else:
        explainer = shap.KernelExplainer(model['model'].predict, shap.kmeans(X_processed, 50))
    shap_values = explainer.shap_values(X_processed)

    # shap plot
    plot_from_shap_values(shap_values, X_processed, save_fig=True, names_list=chem_names_list,
                          plot_title="Metabolite level features",
                          name_file=out_shap_plot, path=OUTPUT_DIR)

    # feature importance
    top_features = feature_importances_from_shap_values(shap_values, X_processed, names_list=chem_names_list, n=20)
top_features.to_csv(OUTPUT_DIR + out_ft_file, index=False)


print("end")
