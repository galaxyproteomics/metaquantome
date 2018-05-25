import pandas as pd
from pandas.plotting import scatter_matrix
import matplotlib.pyplot as plt
import sklearn
from sklearn import model_selection
from sklearn.metrics import classification_report
from sklearn.metrics import confusion_matrix
from sklearn.metrics import accuracy_score
from sklearn.linear_model import LogisticRegression
from sklearn.tree import DecisionTreeClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.naive_bayes import GaussianNB
from sklearn.svm import SVC
import numpy as np

names = ["cog", "peptide", "int737NS", "int737WS", "int852NS", "int852WS", "int867NS", "int867WS"]
df = pd.read_csv("cogs.tabular", sep="\t", index_col="peptide")
df.query("cog == cog", inplace=True)

cog = list(df["cog"])
cog2 = [i.split(',')[0] for i in cog]

df.drop("cog", axis=1, inplace=True)
# df.drop("cog", inplace = True)

log_df = df.apply(np.log, axis = 0)
# scatter_matrix(log_df)
# plt.show()

# create a validation dataset
X = sklearn.preprocessing.scale(df.values)
Y = cog2
validation_size = 0.2
seed = 101
X_train, X_validation, Y_train, Y_validation = model_selection.train_test_split(X, Y, test_size=validation_size, random_state=seed)
scoring = "accuracy"

