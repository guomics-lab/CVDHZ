import numpy as np
from numpy import sort
import pandas as pd
import pdb
from sklearn.model_selection import train_test_split
from scipy.stats import wasserstein_distance
from sklearn.metrics import accuracy_score, roc_auc_score, confusion_matrix
from sklearn.feature_selection import SelectFromModel
from sklearn.ensemble import RandomForestClassifier
from xgboost import XGBClassifier
from sklearn.impute import SimpleImputer
from sklearn.utils import resample
import os
import pickle
import random
import argparse
import csv
import time
from sklearn.model_selection import KFold # import KFold

def cross_validation_split_index_list(n,k,s):
    random.seed(s)
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





def find_best_param_list(X_train,y_train,X_test,y_test,data_split,nindex):
    X_train_1 = X_train[data_split[nindex][0]]
    X_val_1 = X_train[data_split[nindex][1]]
    y_train_1 = y_train[data_split[nindex][0]]
    y_val_1 = y_train[data_split[nindex][1]]

    lr = [0.25,0.255,0.26,0.265,0.27,0.275,0.28,0.285,0.29,0.295,0.3]
    ss = [0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95]
    cbt = [0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95]
    cbl = [0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95]
    spw = [0.2,0.6,1]

    best_measure = 0
    for learning_rate in lr:
        for subsample in ss:
            for colsample_bytree in cbt:
                for colsample_bylevel in cbl:
                    for scale_pos_weight in spw:
                        param_list = [learning_rate,subsample,colsample_bytree,colsample_bylevel,scale_pos_weight]
                        print(param_list)
                        model = XGBClassifier(learning_rate=learning_rate,subsample=subsample,colsample_bytree=colsample_bytree,colsample_bylevel=colsample_bylevel,n_estimators=20,objective='binary:logistic',scale_pos_weight=scale_pos_weight,seed=0,use_label_encoder =False, eval_metric='logloss')

                        model.fit(X_train_1, y_train_1)
                        y_pred = model.predict_proba(X_val_1)[:,1]
                        predictions = [round(value) for value in y_pred]
                        accuracy = accuracy_score(y_val_1, predictions)

                        thresholds = sort(model.feature_importances_)[::-1]
                        count = 4
                        best_all = [0,0,0,0]
                        for thresh in thresholds[4:10]:
                            selection = SelectFromModel(model, threshold=thresh, prefit=True)
                            select_X_train = selection.transform(X_train_1)
                            selection_model = XGBClassifier(n_estimators=20,use_label_encoder =False,
                                eval_metric='logloss')
                            selection_model.fit(select_X_train, y_train_1)
                            select_X_val = selection.transform(X_val_1)
                            y_pred = selection_model.predict_proba(select_X_val)[:,1]
                            predictions = (y_pred > 123 / 137) * 1
                            val_accuracy = accuracy_score(y_val_1, predictions)
                            validation_auc = roc_auc_score(y_val_1, y_pred)
                            count += 1
                            final_data=selection.transform(X_test)
                            y_pred_test = selection_model.predict_proba(final_data)[:,1]
                            predictions_test = (y_pred_test > 123 / 137) * 1
                            #pdb.set_trace()
                            test_accuracy = accuracy_score(y_test, predictions_test)
                            test_auc = roc_auc_score(y_test, y_pred_test)
                            cm = confusion_matrix(y_test, predictions_test, labels=[1, 0])
                            eps = 0.001
                            test_spe = cm[1, 1] / (cm[1, 1] + cm[1, 0] + eps)
                            test_sen = cm[0, 0] / (cm[0, 0] + cm[0, 1] + eps)


                            if validation_auc > 0.999 and validation_auc + val_accuracy > best_measure:
                                best_measure = validation_auc + val_accuracy
                                best_acc = val_accuracy
                                best_auc = validation_auc
                                param_list.append(count)
                                best_param_list = param_list
                                print("Excellent!")
                                print("Val Acc: %.2f, Val AUC: %.2f, Test Acc: %.2f, Test AUC: %.2f, Test SEN: %.2f, Test SPE: %.2f, count: %d" % (val_accuracy * 100.0, validation_auc * 100.0, test_accuracy * 100.0, test_auc * 100.0, test_sen * 100, test_spe * 100, count))
    return best_param_list


if __name__ == "__main__":
    Feature = pd.read_csv('selected_Serum_protein.csv', index_col=0)
    Label = pd.read_csv('labels_all.csv', index_col=0)

    X_train = Feature.loc[Label['group'] != 2].values
    X_test = Feature.loc[Label['group'] == 2].values

    y_train = Label.loc[Label['group'] != 2]['Seroconversion of Nab j , (57d)'].values
    y_test = Label.loc[Label['group'] == 2]['Seroconversion of Nab j , (57d)'].values

    sss = cross_validation_split_index_list(len(y_train), 5, 176) #176

    best_param_list = find_best_param_list(X_train,y_train,X_test,y_test,sss,4)
    print(best_param_list)