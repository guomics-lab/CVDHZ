import numpy as np
from numpy import sort
import pandas as pd
import pdb
from sklearn.model_selection import train_test_split
from scipy.stats import wasserstein_distance
from sklearn.metrics import accuracy_score, roc_auc_score, roc_curve, confusion_matrix
from sklearn.feature_selection import SelectFromModel
from sklearn.ensemble import RandomForestClassifier
from xgboost import XGBClassifier
from sklearn.impute import SimpleImputer
import os
import pickle
import random
from sklearn.model_selection import KFold # import KFold
import heapq
import matplotlib.pyplot as plt
import shap
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})


Feature1 = pd.read_csv('selected_Serum_protein.csv', index_col=0)
Feature2 = pd.read_csv('selected_PBMC_protein.csv', index_col=0)
Label = pd.read_csv('labels_all.csv', index_col=0)

print( Feature1.loc[Label['group'] != 2][['Q9UK55', 'P08294', 'Q13103', 'P54802', 'P09960']].describe())
print( Feature1.loc[Label['group'] == 2][['Q9UK55', 'P08294', 'Q13103', 'P54802', 'P09960']].describe())

X_train1 = Feature1.loc[(Label['group'] != 2) & (Label['Seroconversion of Nab m (28d)'] == 0)]
X_test1 = Feature1.loc[(Label['group'] == 2) & (Label['Seroconversion of Nab m (28d)'] == 0)]

X_train2 = Feature2.loc[(Label['group'] != 2) & (Label['Seroconversion of Nab m (28d)'] == 0)]
X_test2 = Feature2.loc[(Label['group'] == 2) & (Label['Seroconversion of Nab m (28d)'] == 0)]

y_train = Label.loc[(Label['group'] != 2) & (Label['Seroconversion of Nab m (28d)'] == 0)]['Seroconversion of Nab j , (57d)']
y_test = Label.loc[(Label['group'] == 2) & (Label['Seroconversion of Nab m (28d)'] == 0)]['Seroconversion of Nab j , (57d)']

# ------serum data input--------
serum_selected_feats = ['Q969E1', 'P49747', 'Q13103', 'Q9UK55', 'P13497']
select_X_train1 = Feature1.loc[(Label['group'] != 2) & (Label['Seroconversion of Nab m (28d)'] == 0)][serum_selected_feats]
select_X_test1 = Feature1.loc[(Label['group'] == 2) & (Label['Seroconversion of Nab m (28d)'] == 0)][serum_selected_feats]

selection_model1 = XGBClassifier(colsample_bylevel=0.5, learning_rate=0.1, max_bin=5, max_depth=2, min_child_weight=3, \
    n_estimators=30, random_state=1000, scale_pos_weight=0.3, subsample=0.8, gamma=0)
selection_model1.fit(select_X_train1, y_train)

explainer1 = shap.TreeExplainer(selection_model1)
shap_values1 = explainer1.shap_values(select_X_train1)
data_for_shap1 = pd.DataFrame(select_X_train1, columns=serum_selected_feats)
shap.summary_plot(shap_values1, data_for_shap1, show=False)

# shap.summary_plot(shap_values1, data_for_shap1, plot_type="bar")
plt.savefig('shap_Serum.pdf')
plt.cla()

y_pred_test1 = selection_model1.predict_proba(select_X_test1)[:,1]
predictions_test1 = [round(value) for value in y_pred_test1]
val_auc1 = roc_auc_score(y_test, y_pred_test1)
print('---Serum TEST---', val_auc1)

# ------pbmc data input--------
pbmc_select_features = ['Q6P1L8', 'Q9H3U1', 'Q9NUJ1', 'O00443', 'P01871']
select_X_train2 = X_train2[pbmc_select_features]
select_X_test2 = X_test2[pbmc_select_features]
selection_model2 = XGBClassifier(colsample_bylevel=0.3, learning_rate=0.05, max_bin=5, max_depth=3, min_child_weight=3, \
    n_estimators=200, random_state=1000, scale_pos_weight=1.0, subsample=1.0, gamma=0, colsample_bytree=0.6)
selection_model2.fit(select_X_train2, y_train)

explainer2 = shap.TreeExplainer(selection_model2)
shap_values2 = explainer2.shap_values(select_X_train2)
data_for_shap2 = pd.DataFrame(select_X_train2, columns=pbmc_select_features)
shap.summary_plot(shap_values2, data_for_shap2, show=False)
# shap.summary_plot(shap_values2, data_for_shap2, plot_type="bar")
plt.savefig('shap_PBMC.pdf')
plt.cla()

y_pred_test2 = selection_model2.predict_proba(select_X_test2)[:,1]
predictions_test2 = [round(value) for value in y_pred_test2]
val_auc2 = roc_auc_score(y_test, y_pred_test2)
print('---PMBC TEST', val_auc2)

all_test_y_preds = (y_pred_test1 + y_pred_test2)/2
all_test_predictions = (all_test_y_preds>0.6)*1
test_accuracy = accuracy_score(y_test, all_test_predictions)

test_auc = roc_auc_score(y_test, all_test_y_preds)
print("Test Accuracy: %.2f%%" % (test_accuracy*100.0))
print("Test AUC: %.2f%%" % (test_auc*100.0))

cm = confusion_matrix(y_test,all_test_predictions,labels=[1, 0])
print("TN is %d" % cm[0,0])
print("FP is %d" % cm[0,1])
print("FN is %d" % cm[1,0])
print("TP is %d" % cm[1,1])

file_pred = pd.DataFrame(index=Label.loc[(Label['group'] == 2) & (Label['Seroconversion of Nab m (28d)'] == 0)]['Seroconversion of Nab j , (57d)'].index,columns=['label','PBMC','Serum','combine'])
file_pred['label'] = y_test
file_pred['PBMC'] = y_pred_test2
file_pred['Serum'] = y_pred_test1
file_pred['combine'] = all_test_y_preds
file_pred.to_csv('pred_score.csv')

fpr_val1,tpr_val1,threshold_val1 = roc_curve(y_test, y_pred_test1)
fpr_val2,tpr_val2,threshold_val2 = roc_curve(y_test, y_pred_test2)
fpr_test,tpr_test,threshold_test = roc_curve(y_test, all_test_y_preds)
lw = 2
plt.plot(fpr_val1, tpr_val1, color='darkorange',lw=lw, label='Serum ROC curve (area = %0.2f)' % val_auc1)
plt.plot(fpr_val2, tpr_val2, color='blue',lw=lw, label='PBMC ROC curve (area = %0.2f)' % val_auc2)
plt.plot(fpr_test, tpr_test, color='red',lw=lw, label='Average ROC curve (area = %0.2f)' % test_auc)
plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Receiver operating characteristic example')
plt.legend(loc="lower right")
plt.savefig(r'ALL_ROC_multi.pdf')
plt.cla()

colors = ['b','r']
for i in range(len(all_test_y_preds)):
    if (all_test_y_preds[i]<0.6 and y_test[i]< 0.6) or (all_test_y_preds[i]>=0.6 and y_test[i]< 0.6):
        k = 0
    else:
        k = 1
    plt.scatter(all_test_y_preds[i], i, c=colors[k])
plt.vlines(123/137, 0, len(all_test_y_preds), colors = "black", linestyles = "dashed")
plt.title('ALL_scatter_multi')
plt.savefig(r'ALL_scatter_multi.pdf')
plt.cla()

colors = ['b','r']
for i in range(len(y_pred_test2)):
    if (y_pred_test2[i]<0.5 and y_test[i]< 0.5) or (y_pred_test2[i]>=0.5 and y_test[i]< 0.5):
        k = 0
    else:
        k = 1
    plt.scatter(y_pred_test2[i], i, c=colors[k])
plt.vlines(123/137, 0, len(y_pred_test2), colors = "black", linestyles = "dashed")
plt.title('PBMC')
plt.savefig(r'ALL_scatter_PBMC.pdf')
plt.cla()

colors = ['b','r']
for i in range(len(y_pred_test1)):
    if (y_pred_test1[i]<0.5 and y_test[i]< 0.5) or (y_pred_test1[i]>=0.5 and y_test[i]< 0.5):
        k = 0
    else:
        k = 1
    plt.scatter(y_pred_test1[i], i, c=colors[k])
plt.vlines(123/137, 0, len(y_pred_test1), colors = "black", linestyles = "dashed")
plt.title('Serum')
plt.savefig(r'ALL_scatter_Serum.pdf')
plt.cla()

import sys
sys.exit(0)
feature_set = set(max_feature_index_list1)
feature_set = list(feature_set)
feature_set_score = [0] * len(feature_set)
sub_feature_set = [max_feature_index_list1]
score_list = [p1]

for ii in range(len(feature_set)):
    for kk in score_list:
        if kk[feature_set[ii]] > feature_set_score[ii]:
            feature_set_score[ii] = kk[feature_set[ii]]
score_rank = np.argsort(feature_set_score)[::-1]
feature_rank = np.array(feature_set)[score_rank]
feature_score_rank = np.array(feature_set_score)[score_rank]
X_all_colname = pd.read_csv('selected_PBMC_protein.csv', index_col=0).columns
extracted_colname = X_all_colname[feature_rank]
plt.barh(range(len(extracted_colname)), feature_score_rank[::-1],color='b',tick_label=extracted_colname[::-1])
plt.tick_params(labelsize=6)
plt.tight_layout()
plt.savefig(r'PBMC_feature_importance.pdf')
plt.cla()
