"""
Created on Thu Feb 01 14:19:06 2024
Updated on Sun Apr 07 18:14:32 2024
@author: vkeggers and virallyDanny
@description:

@Usage: `python xgboost.py --Genus <Genus>`
"""

import pandas as pd
import numpy as np
import xgboost
from xgboost import XGBClassifier
from sklearn.model_selection import train_test_split, KFold, cross_val_score
from sklearn.metrics import classification_report
import matplotlib.pyplot as plt
import argparse
import warnings
from sklearn.metrics import roc_curve, auc
warnings.filterwarnings('ignore')

# Create an argument parser object to parse the value of the --Genus argument to the genus variable
# parser = argparse.ArgumentParser(description='XGBoost model to predict the origin of a TE')
# parser.add_argument('--Genus', type=str, help='Genus of the TE')
# args = parser.parse_args()
# genus = args.Genus

genus = 'Oscheius'

# Load SIDR stats tsv file to a pandas dataframe titled data
stats = pd.read_csv('./data/SIDRstats.tsv', sep='\t', header=0)

# Make a new dataframe titled trainingDF that only contains rows which do not equal 'No hits found' in the 'Origin' column
trainingDF = stats[stats['Origin'] != 'No hits found']

# Make a new dataframe titled orig where if the 'Origin' column equals the 'Genus' string variable. If it is True then change the Origin value to True, else False
orig = trainingDF[['Origin']]
orig['Origin'] = trainingDF['Origin'].str.contains(genus).astype(int)

# Make a variable titled "Train" that samples trainingDF with a random sampling 1/3 of the data in the dataframe
Train = trainingDF.sample(frac=1/2)

# Make a variable titled stat1_test which takes all the rows not in Train and keeps all columns except the 'contig' and 'Origin' columns
stat1_test = trainingDF[~trainingDF.index.isin(Train.index)].drop(['contig', 'Origin'], axis=1)

stat1_train = trainingDF[trainingDF.index.isin(Train.index)].drop(['contig', 'Origin'], axis=1)

origin_test = orig[~orig.index.isin(Train.index)]

origin_train = orig[orig.index.isin(Train.index)]

####################
### WORKS SO FAR ###
####################

model = XGBClassifier(n_estimator=1000, max_depth=6, reg_lambda=2, random_state=3, objective='binary:logistic', distribution='bernoulli')

model.fit(stat1_train, origin_train)

kfold = KFold(n_splits=5, random_state=7, shuffle=True)

results = cross_val_score(model, stat1_train, origin_train, cv=kfold)

print("Accuracy: %.2f%% (%.2f%%)" % (results.mean()*100, results.std()*100))
model.save_model("model-1.json")

# Print SIDR predictions to stats dataframe
stats['SIDR_predictions'] = model.predict(stats.drop(['contig', 'Origin'], axis=1))

# Save the stats dataframe to a new tsv file
stats.to_csv('./results/Full_SIDR_predictions.tsv', sep='\t', index=False)

# Save the stats dataframe to a new tsv file if the SIDR_predictions column equals 1
stats[stats['SIDR_predictions'] == 1].to_csv('./results/Kept_SIDR_predictions.tsv', sep='\t', index=False)
stats[stats['SIDR_predictions'] == 0].to_csv('./results/Removed_SIDR_predictions.tsv', sep='\t', index=False)

# Plot the decision tree
ax = xgboost.plot_tree(model, num_trees=1)
plt.savefig('./figures/tree_plot.png', dpi=1200, bbox_inches='tight')
plt.close()

# Feature Importance
feature_importance = model.feature_importances_features = stats.columns

# Plot the top 7 features
xgboost.plot_importance(model, max_num_features=12)

# Prediction Report
y_pred = model.predict(stat1_test)
report = classification_report(origin_test, y_pred)

# Show the plot
plt.show()
plt.savefig('./figures/feature_importance.png', dpi=1200, bbox_inches='tight')
plt.close()


# Get predicted probabilities for positive class
y_prob = model.predict_proba(stat1_test)[:, 1]

# Compute ROC curve and ROC area
fpr, tpr, _ = roc_curve(origin_test, y_prob)
roc_auc = auc(fpr, tpr)

# Plot ROC curve
plt.figure()
lw = 2
plt.plot(fpr, tpr, color='darkorange', lw=lw, label='ROC curve (area = %0.2f)' % roc_auc)
plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Receiver Operating Characteristic Curve')
plt.legend(loc="lower right")
plt.show()
plt.savefig('./figures/ROCcurve.png', dpi=1200, bbox_inches='tight')


# Get predicted probabilities for positive class for training and test sets
train_y_prob = model.predict_proba(stat1_train)[:, 1]
test_y_prob = model.predict_proba(stat1_test)[:, 1]

# Compute ROC curve and ROC area for training set
train_fpr, train_tpr, _ = roc_curve(origin_train, train_y_prob)
train_roc_auc = auc(train_fpr, train_tpr)

# Compute ROC curve and ROC area for test set
test_fpr, test_tpr, _ = roc_curve(origin_test, test_y_prob)
test_roc_auc = auc(test_fpr, test_tpr)

# Plot ROC curve for training set
plt.figure()
lw = 2
plt.plot(train_fpr, train_tpr, color='darkorange', lw=lw, label='Training ROC curve (area = %0.2f)' % train_roc_auc)

# Plot ROC curve for test set
plt.plot(test_fpr, test_tpr, color='blue', lw=lw, label='Test ROC curve (area = %0.2f)' % test_roc_auc)

plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Receiver Operating Characteristic Curve')
plt.legend(loc="lower right")
plt.show()
plt.savefig('./figures/JointROCcurve.png', dpi=1200, bbox_inches='tight')
