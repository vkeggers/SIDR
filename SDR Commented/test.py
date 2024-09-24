# Import necessary libraries from scikit-learn for the Decision Tree Classifier and model selection utilities.
from sklearn import tree
from sklearn.tree import DecisionTreeClassifier
from sklearn import model_selection as ms

# Initialize lists to store sequence data categorized as target, not target, and unclassified.
target = []
not_target = []
unclassified = []

# Function to import sequence data and categorize it based on certain criteria.
def import_seqdata(seqname, seqlen, gc, cov, hasBlast, isTarget, kdist_list):
    # Create a tuple for the sequence data with the specified parameters.
    seqdata = (seqname, seqlen, kdist_list, isTarget, gc, cov)
    # Categorize the sequence data based on the presence of Blast hits and target status.
    if hasBlast == 1:
        if isTarget == 1:
            target.append(seqdata)
        else:
            not_target.append(seqdata)
    else:
        unclassified.append(seqdata)
    return 0

# Function to classify unclassified sequences using a pre-trained Decision Tree classifier.
def DT_classify(classifier):
    regionIDs = []  # To store sequence names.
    X = []  # To store features for prediction.
    # Extract relevant features from unclassified sequences.
    for item in unclassified:
        regionIDs.append(item[0])  # region name
        X.append(item[4:])
    Y = []  # To store prediction results.
    # Use the classifier to predict the category of each unclassified sequence.
    for i in classifier.predict(X):
        Y.append(i)
    # Pair each region ID with its predicted category.
    return list(zip(regionIDs, Y))

# Function to build and train the Decision Tree model, then classify unclassified sequences.
def DT_model(unclassified):
    X = []  # Feature set.
    Y = []  # Target variable (isTarget).
    # Combine target and not target data for training.
    for item in target + not_target:
        X.append(item[4:])  # GC content and coverage.
        Y.append(item[3])  # isTarget.
    # Split data into training and testing sets.
    X_train, X_test, Y_train, Y_test = ms.train_test_split(X, Y, test_size=0.33, random_state=0)
    # Initialize and train the Decision Tree classifier.
    classifier = DecisionTreeClassifier().fit(X_train, Y_train)
    # Evaluate and print the classifier's performance on the test set.
    print("Classifier built, score is %s out of 1.00" % classifier.score(X_test, Y_test))
    # Return the result of classifying unclassified sequences with the trained classifier.
    return DT_classify(classifier)

# Function to run the overall analysis, import data, train the model, and classify sequences.
def run_analysis(kmer_list):
    print("Seqdata imported, %d target sequences, %d contaminant sequences, and %d unclassified sequences\n" % (len(target), len(not_target), len(unclassified)))
    # Proceed with model training and classification if there are unclassified sequences.
    if len(unclassified) != 0:
        result = DT_model(unclassified)
    else:
        result = []
    # Write classification results to separate files for further analysis or review.
    with open("tokeep.txt", "w+") as k, open("toremove.txt", "w+") as r:
        for item in result:
            if item[1] == 1:  # Target sequences to keep.
                k.write("%s x\n" % item[0])
            else:  # Non-target sequences to remove.
                r.write("%s x\n" % item[0])
    # Optionally, write categorized target and non-target sequences to files.
    with open("target.txt", "w+") as t, open("nottarget.txt", "w+") as n:
        for item in target:
            t.write("%s x\n" % item[0])
        for item in not_target:
            n.write("%s x\n" % item[0])
    return 1
 in target:
		t.write("%s x\n" % item)
	for item in nottarget:
		n.write("%s x\n" % item)


    return 1
