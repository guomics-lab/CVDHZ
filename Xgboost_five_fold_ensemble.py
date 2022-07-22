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


def cross_validation_split_index_list(n,k):
    random.seed(176)
    index = np.arange(n)
    random.shuffle(index)
    list = []
    first_index = 0
    step = round(n/k)
    for i in range(k):
        sub_list = []
        if i != k-1:
            test_temp_index = index[first_index:(first_index+step)]
            diff = set(index).difference(set(test_temp_index))
            train_temp_index = []
            for j in diff:
                train_temp_index.append(j)
        else:
            test_temp_index = index[first_index:]
            diff = set(index).difference(set(test_temp_index))
            train_temp_index = []
            for j in diff:
                train_temp_index.append(j)
        sub_list.append(train_temp_index)
        sub_list.append(test_temp_index)
        list.append(sub_list)
        first_index += step
    return list


Feature1 = pd.read_csv('selected_Serum_protein.csv', index_col=0)
Feature2 = pd.read_csv('selected_PBMC_protein.csv', index_col=0)
Label = pd.read_csv('labels_all.csv', index_col=0)

X_train1 = Feature1.loc[Label['group'] != 2].values
X_test1 = Feature1.loc[Label['group'] == 2].values

X_train2 = Feature2.loc[Label['group'] != 2].values
X_test2 = Feature2.loc[Label['group'] == 2].values


y_train = Label.loc[Label['group'] != 2]['Seroconversion of Nab j , (57d)'].values
y_test = Label.loc[Label['group'] == 2]['Seroconversion of Nab j , (57d)'].values
sss = cross_validation_split_index_list(len(y_train),5)
X_train_1 = X_train1[sss[4][0]]
X_val_1 = X_train1[sss[4][1]]
y_train_1 = y_train[sss[4][0]]
y_val_1 = y_train[sss[4][1]]
model1 = XGBClassifier(learning_rate=0.265,subsample=0.55,colsample_bytree=0.6,colsample_bylevel=0.95,n_estimators=20,objective='binary:logistic',scale_pos_weight=0.2,seed=0,use_label_encoder =False, eval_metric='logloss')
model1.fit(X_train_1, y_train_1)
thresh1 = sort(model1.feature_importances_)[::-1][6]
p1 = list(model1.feature_importances_)
max_feature_index_list1 = [*map(p1.index, heapq.nlargest(7, p1))]
selection1 = SelectFromModel(model1, threshold=thresh1, prefit=True)
select_X_train1 = selection1.transform(X_train_1)
selection_model1 = XGBClassifier(n_estimators=20,use_label_encoder =False, eval_metric='logloss',objective='binary:logistic')
selection_model1.fit(select_X_train1, y_train_1)

explainer1 = shap.TreeExplainer(selection_model1)
shap_values1 = explainer1.shap_values(select_X_train1)
col_for_shap1 = selection1.transform(Feature1.loc[Label['group'] != 2].columns.values.reshape(1,-1))[0]
data_for_shap1 = pd.DataFrame(select_X_train1,columns=col_for_shap1)
shap.summary_plot(shap_values1, data_for_shap1, plot_type="bar")
plt.savefig('shap_Serum.pdf')
plt.cla()

select_X_val1 = selection1.transform(X_val_1)
y_pred1 = selection_model1.predict_proba(select_X_val1)[:,1]
predictions1 = [round(value) for value in y_pred1]
select_X_test1=selection1.transform(X_test1)
y_pred_test1 = selection_model1.predict_proba(select_X_test1)[:,1]
predictions_test1 = [round(value) for value in y_pred_test1]
val_auc1 = roc_auc_score(y_test, y_pred_test1)
print(val_auc1)

X_train_2 = X_train2[sss[2][0]]
X_val_2 = X_train2[sss[2][1]]
y_train_2 = y_train[sss[2][0]]
y_val_2 = y_train[sss[2][1]]
pdb.set_trace()
model2 = XGBClassifier(learning_rate=0.29,subsample=0.8,colsample_bytree=0.95,colsample_bylevel=0.6,n_estimators=20,objective='binary:logistic',scale_pos_weight=1,seed=0,use_label_encoder =False, eval_metric='logloss')
model2.fit(X_train_2, y_train_2)
explainer2 = shap.TreeExplainer(model2)
thresh2 = sort(model2.feature_importances_)[::-1][4]
p2 = list(model2.feature_importances_)
max_feature_index_list2 = [*map(p2.index, heapq.nlargest(5, p2))]
selection2 = SelectFromModel(model2, threshold=thresh2, prefit=True)
select_X_train2 = selection2.transform(X_train_2)
selection_model2 = XGBClassifier(n_estimators=20,use_label_encoder =False, eval_metric='logloss',objective='binary:logistic')
selection_model2.fit(select_X_train2, y_train_2)

explainer2 = shap.TreeExplainer(selection_model2)
shap_values2 = explainer2.shap_values(select_X_train2)
col_for_shap2 = selection2.transform(Feature2.loc[Label['group'] != 2].columns.values.reshape(1,-1))[0]
data_for_shap2 = pd.DataFrame(select_X_train2,columns=col_for_shap2)
shap2 = shap.summary_plot(shap_values2, data_for_shap2, show=True)
shap.summary_plot(shap_values2, data_for_shap2, plot_type="bar")
plt.savefig('shap_PBMC.pdf')
plt.cla()

select_X_val2 = selection2.transform(X_val_2)
y_pred2 = selection_model2.predict_proba(select_X_val2)[:,1]
predictions2 = [round(value) for value in y_pred2]
select_X_test2=selection2.transform(X_test2)
y_pred_test2 = selection_model2.predict_proba(select_X_test2)[:,1]
predictions_test2 = [round(value) for value in y_pred_test2]
val_auc2 = roc_auc_score(y_test, y_pred_test2)
print(val_auc2)

# val_auc = roc_auc_score(y_val_2, y_pred2)
# test_auc = roc_auc_score(y_test, y_pred_test2)
# fpr_val,tpr_val,threshold_val = roc_curve(y_val_2, y_pred2)
# fpr_test,tpr_test,threshold_test = roc_curve(y_test, y_pred_test2)
# lw = 2
# plt.plot(fpr_val, tpr_val, color='darkorange',lw=lw, label='Validation ROC curve (area = %0.2f)' % val_auc)
# plt.plot(fpr_test, tpr_test, color='red',lw=lw, label='Test ROC curve (area = %0.2f)' % test_auc)
# plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
# plt.xlim([0.0, 1.0])
# plt.ylim([0.0, 1.05])
# plt.xlabel('False Positive Rate')
# plt.ylabel('True Positive Rate')
# plt.title('Receiver operating characteristic example')
# plt.legend(loc="lower right")
# plt.savefig(r'ALL_ROC_multi.pdf')
# plt.cla()

all_test_y_preds = (y_pred_test1+y_pred_test2)/2
all_test_predictions = (all_test_y_preds>123/137)*1
test_accuracy = accuracy_score(y_test, all_test_predictions)

test_auc = roc_auc_score(y_test, all_test_y_preds)
print("Test Accuracy: %.2f%%" % (test_accuracy*100.0))
print("Test AUC: %.2f%%" % (test_auc*100.0))

cm = confusion_matrix(y_test,all_test_predictions,labels=[1, 0])
print("TN is %d" % cm[0,0])
print("FP is %d" % cm[0,1])
print("FN is %d" % cm[1,0])
print("TP is %d" % cm[1,1])

file_pred = pd.DataFrame(index=Label.loc[Label['group'] == 2]['Seroconversion of Nab j , (57d)'].index,columns=['label','PBMC','Serum','combine'])
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
    if (all_test_y_preds[i]<0.5 and y_test[i]< 0.5) or (all_test_y_preds[i]>=0.5 and y_test[i]< 0.5):
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
