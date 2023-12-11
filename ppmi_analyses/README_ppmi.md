Scripts utilized for ML modelling of transcritpomics data from the PPMI cohort for mild vs. severe UPDRS III sample classification.



### Data pre-processing, prior to ML modelling

##### ppmi_generate_temp_features.R
Generate longitudinal features (lm-time, sd) from longitudinal transcriptomcis data (T0, T1, T2, T3).

##### ppmi_extract_data_updrs3BL.R, ppmi_extract_data_updrs3V08.R
Filters phenotypic and transcriptommics data of PD patients having UPDRS III score at baseline (BL), and PD patients with UPDRS III at the fourth timepoint (V08) along with extracted temporal features, respectively. 

##### ppmi_makebinaryvar.R
Transform UPDRS__3 continuous variable to binary (0,1) correspondingly to (below median, above median) classes.

##### ppmi_data4ML_class.R, ppmi_data4ML_TS_class.R
Perform unsupervised filters to generate data for ML modelling of mild/severe UPDRS III from RNAseq snapshot data (T0), and longitudinal features extracted from longitudinal RNAseq data, respectively.

* These scripts employ as input transcriptomics and phenotypical data resulting from previous pre-processing scripts described in repository [statistical_analyses_cross_long_PD](https://gitlab.lcsb.uni.lu/elisa.gomezdelope/statistical_analyses_cross_long_pd) for **parsing data** and **Baseline (T0) PD/HC** (ppmi_filter_gene_expression.R, ppmi_norm_gene_expression.R, ppmi_generate_pathway_level.R, ppmi_deseq_salmon_star.R). 



### ML modelling with LASSO feature selection at baseline and using temporal features

Please, install the environment package *digipd_ml* with the class **NestedCV()** to run these scripts.

##### PD_control_train.py, PD_control_train_PW.py
Perform training for mild/severe UPDRS III classification with multiple binary classifiers on transcriptomics data in a nested cross-validation setting, allowing for undersampling, feature scaling, and feature selection with LASSO. Input taken from *ppmi_data4ML_class.R*, *ppmi_data4ML_TS_class.R* outputs (which include the unsupervised feature selection filters).

##### PD_control_finalmodel_plot.py, PD_control_finalmodel_plot_PW.py
Train a selected model on the entire training set for mild/severe UPDRS III classification and validate on the test set, compute SHAP values and generate corresponding plots and tables.

PD_control_finalmodel_plot_notest.R, PD_control_finalmodel_plot_notest_PW.py
Train a selected model on the entire training set for mild/severe UPDRS III classification, compute SHAP values and generate corresponding plots and tables for such model.

##### DB_st_plots.ipynb, plots.py
Generate figures from the nested CV results on different dynamic features and baseline.



### Post-hoc analysis

##### ppmi_friedman_UPDRS3.R
Performs comparisons among model performance, data types performance, and feature selection methods using the friedman test with Bergmann-Hommel corrections for pairwise comparisons, and with Holm corrections for one-to-many comparisons.



