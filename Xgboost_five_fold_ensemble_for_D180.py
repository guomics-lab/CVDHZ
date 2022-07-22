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
    random.seed(202)
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


Feature = pd.read_csv('all_feature_PBMC.csv', index_col=0)
Label = pd.read_csv('clinicinfo_127samples_matched.csv', index_col=0)

X_train = Feature.loc[Label.loc[Label['Recruitment_date'] != 2].index].values[:,:52]
X_test = Feature.loc[Label.loc[Label['Recruitment_date'] == 2].index].values[:,:52]
#X_test2 = Feature.loc[Label['group'] == 2].values


#pdb.set_trace()
y_train = Label.loc[Label['Recruitment_date'] != 2]['180d-anti'].values
y_test = Label.loc[Label['Recruitment_date'] == 2]['180d-anti'].values

sss = cross_validation_split_index_list(len(y_train),5)
#pdb.set_trace()
X_train_1 = X_train[sss[0][0]]
X_val_1 = X_train[sss[0][1]]
y_train_1 = y_train[sss[0][0]]
y_val_1 = y_train[sss[0][1]]
model1 = XGBClassifier(learning_rate=0.3,subsample=0.6,colsample_bytree=0.7,colsample_bylevel=0.5,n_estimators=20,objective='binary:logistic',scale_pos_weight=0.6,seed=9,use_label_encoder =False, eval_metric='logloss')
model1.fit(X_train_1, y_train_1)
thresh1 = sort(model1.feature_importances_)[::-1][3]
p1 = list(model1.feature_importances_)
max_feature_index_list1 = [*map(p1.index, heapq.nlargest(4, p1))]
selection1 = SelectFromModel(model1, threshold=thresh1, prefit=True)
select_X_train1 = selection1.transform(X_train_1)
selection_model1 = XGBClassifier(n_estimators=20,use_label_encoder =False, eval_metric='logloss',seed=3,objective='binary:logistic')
selection_model1.fit(select_X_train1, y_train_1)

explainer1 = shap.TreeExplainer(selection_model1)
shap_values1 = explainer1.shap_values(select_X_train1)
#pdb.set_trace()
feature_names = selection1.transform(Feature.columns.values[:52].reshape(1,-1))[0]
data_for_shap1 = pd.DataFrame(select_X_train1,columns=feature_names)
#shap1 = shap.summary_plot(shap_values1, data_for_shap1, show=False)
shap.summary_plot(shap_values1, data_for_shap1, plot_type="bar")
# plt.savefig('shap_values.pdf')
# pdb.set_trace()

select_X_val1 = selection1.transform(X_val_1)
y_pred1 = selection_model1.predict_proba(select_X_val1)[:,1]
predictions1 = [round(value) for value in y_pred1]
select_X_test1=selection1.transform(X_test)
y_pred_test1 = selection_model1.predict_proba(select_X_test1)[:,1]
predictions_test1 = [round(value) for value in y_pred_test1]

val_auc = roc_auc_score(y_val_1, y_pred1)
test_auc = roc_auc_score(y_test, y_pred_test1)
import matplotlib.pyplot as plt
fpr_val,tpr_val,threshold_val = roc_curve(y_val_1, y_pred1)
fpr_test,tpr_test,threshold_test = roc_curve(y_test, y_pred_test1)
lw = 2
plt.plot(fpr_val, tpr_val, color='darkorange',lw=lw, label='Validation ROC curve (area = %0.2f)' % val_auc)
plt.plot(fpr_test, tpr_test, color='red',lw=lw, label='Test ROC curve (area = %0.2f)' % test_auc)
plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Receiver operating characteristic example')
plt.legend(loc="lower right")
#plt.savefig(r'ALL_ROC_multi.pdf')
plt.show()


