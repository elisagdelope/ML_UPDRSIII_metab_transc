Scripts utilized for ML modelling of metabolomics data from the luxPARK cohort for mild vs. severe UPDRS III sample classification.


### Data pre-processing, prior to ML modelling

##### lx_generate_temp_features.R
Generate longitudinal features (lm-time, sd) from longitudinal metabolomics data (T0, T1, T2, T3).

##### lx_extract_data_updrs3V0.R, lx_extract_data_updrs3V2.R
Filter phenotypic and metabolomics data of PD patients having UPDRS III score at baseline (V0), and PD patients with UPDRS III at the third timepoint (V2) along with extracted temporal features, respectively. 

##### lx_makebinaryvar.R
Transform UPDRS__3 continuous variable to binary (0,1) correspondingly to (below median, above median) classes.

##### lx_data4ML_class.R, lx_data4ML_TS_class.R
Perform unsupervised filters to generate data for ML modelling of metabolomics snapshot data (T0), and longitudinal features extracted from longitudinal metabolomics data, respectively.

* These scripts employ as input metabolomics and phenotypical data resulting from previous pre-processing scripts described in repository [statistical_analyses_cross_long_PD](https://gitlab.lcsb.uni.lu/elisa.gomezdelope/statistical_analyses_cross_long_pd) for **parsing data** and **Baseline (T0) PD/HC** (lx_extract_visit.R, lx_denovo_filter.R, lx_generate_pathway_level.R). 



### ML modelling with LASSO feature selection at baseline and using temporal features

Please, install the environment package *digipd_ml* with the class **NestedCV()** to run these scripts.

##### PD_control_train.py, PD_control_train_PW.py
Perform training for mild/severe UPDRS III classification with multiple binary classifiers on metabolomics data in a nested cross-validation setting, allowing for undersampling, feature scaling, and feature selection with LASSO. Input taken from *lx_data4ML_class.R*, outputs (which include the unsupervised feature selection filters).

##### PD_control_finalmodel_plot.py, PD_control_finalmodel_plot_PW.py
Train a selected model on the entire training set for mild/severe UPDRS III classification and validate on the test set, compute SHAP values and generate corresponding plots and tables.

PD_control_finalmodel_plot_notest.R, PD_control_finalmodel_plot_notest_PW.py
Train a selected model on the entire training set for mild/severe UPDRS III classification, compute SHAP values and generate corresponding plots and tables for such model.

##### DB_st_plots.ipynb, plots.py
Generate figures from the nested CV results on different dynamic features and baseline.



### Post-hoc analysis

##### lx_friedman_UPDRS3.R
Performs comparisons among model performance, data types performance, and feature selection methods using the friedman test with Bergmann-Hommel corrections for pairwise comparisons, and with Holm corrections for one-to-many comparisons.

