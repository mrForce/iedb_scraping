import matplotlib.pyplot as plt
import argparse
from sklearn.metrics import roc_curve, auc
import numpy as np
import pickle
from TrainingSystems import classification_pwm, classification_hmm
from os import path
from Bio.Alphabet import IUPAC
from itertools import cycle

from sklearn import svm, datasets

from sklearn.model_selection import StratifiedKFold
from enum import Enum
import sys
class PredictionMethod(Enum):
    PWM = 1
    HMM = 2
    QuantMatrix = 3


def compute_mean_roc(fprs, tprs):
    """
     figure out the FPRs where a ROC curve jumps.
    """
    fpr_samples = np.array(sorted(list(set().union(*fprs))))
    #so, sample each FPR/TPR at fpr_samples using interpolation.
    mean_samples = sum([np.interp(fpr_samples, fpr, tpr) for fpr, tpr in zip(fprs, tprs)])/len(fprs)
    """
    Okay, so we've sampled at the points we need to sample at. 
    If we plotted what we have now, we would end up with linear curves between each jump point. We don't want this.

    So, if we have a curve with x values[0.0, 0.1, 0.5, 0.6, 1.0] and y values [0.0, 0.2, 0.7, 0.8, 1.0]
    We need to add the points (0.1, 0.0), (0.5, 0.2), (
    """
    x_values = [fpr_samples[0]]
    y_values = [mean_samples[0]]
    for i in range(1, len(mean_samples)):
        x_values.append(fpr_samples[i])
        x_values.append(fpr_samples[i])
        y_values.append(mean_samples[i - 1])
        y_values.append(mean_samples[i])
    return (x_values, y_values)

"""
data should be of the form (sequence, measurement)

All sequences with measurements less than or equal to  measurement_threshold are binders.
All sequences with measumerents greater than measurement_threshold are non-binders.
"""
def roc_curve_cross(data, measurement_threshold, title, save_location, method, num_rows = False):
    #first, convert the measurements to True and False values based on the measurement threshold
    sequences = np.array([sequence for sequence, measurement in data])
    binary_results = np.array([1 if measurement <= measurement_threshold else 0 for sequence,measurement in data])
    cross_validation = StratifiedKFold(n_splits=5)

    background_sequences = list()
    with open('nine_chunks.txt', 'r') as g:
        background_sequences = list(filter(lambda x: len(x) == 9, [x.strip() for x in g.read().split('\n')]))
    """
    So, for the mean ROC curve, we sample the FPR and TPR at 100 evenly spaced points. That's what they did here: http://scikit-learn.org/stable/auto_examples/model_selection/plot_roc_crossval.html They used interpolation -- I think this is good idea, as it helps to fill in what happened between 

    I was thinking of re-sampling at every FPR where the TPR of a 
    """
    back_fprs = list()
    back_tprs = list()
    no_back_fprs = list()
    no_back_tprs = list()
    j = 1
    for training_indices, test_indices in cross_validation.split(sequences, binary_results):
        if method == PredictionMethod.PWM:
            background_binary = classification_pwm.PWMBinaryClassifier(background_sequences)
            no_background_binary = classification_pwm.PWMBinaryClassifier()
        if method == PredictionMethod.HMM:
            background_binary = classification_hmm.HMMBinaryClassifier(num_rows, background_sequences)
            no_background_binary = classification_hmm.HMMBinaryClassifier(num_rows)
        training_data_sequences = sequences[training_indices]
        training_data_binary_results = binary_results[training_indices]
        training_data = list(zip(training_data_sequences, training_data_binary_results))
        background_binary.train(training_data, IUPAC.protein.letters)
        no_background_binary.train(training_data, IUPAC.protein.letters)
        #we subtract since we are working with log probabilities.
        if method == PredictionMethod.PWM:
            background_test_scores = np.array([background_binary.get_positive_score(peptide) - background_binary.get_negative_score(peptide) for peptide in sequences[test_indices]])
            no_background_test_scores = np.array([no_background_binary.get_positive_score(peptide) - no_background_binary.get_negative_score(peptide) for peptide in sequences[test_indices]])
            #test_scores = np.array([binary_classifier.get_positive_score(peptide) for peptide in sequences[test_indices]])
        if method == PredictionMethod.HMM:
            background_test_scores = np.array([background_binary.get_positive_score(peptide) for peptide in sequences[test_indices]])
            no_background_test_scores = np.array([no_background_binary.get_positive_score(peptide) for peptide in sequences[test_indices]])
        print('test indices')
        print(binary_results[test_indices])


        back_fpr, back_tpr, back_thresholds = roc_curve(binary_results[test_indices], background_test_scores)
        no_back_fpr, no_back_tpr, no_back_thresholds = roc_curve(binary_results[test_indices], no_background_test_scores)
        back_fprs.append(back_fpr)
        back_tprs.append(back_tpr)
        no_back_fprs.append(no_back_fpr)
        no_back_tprs.append(no_back_tpr)
        back_auc = auc(back_fpr, back_tpr)
        no_back_auc = auc(no_back_fpr, no_back_tpr)
        j+= 1
    back_mean_fpr, back_mean_tpr = compute_mean_roc(back_fprs, back_tprs)
    back_mean_area = auc(back_mean_fpr, back_mean_tpr)
    plt.plot(back_mean_fpr, back_mean_tpr, lw=2, color='black', linestyle='--', label='Background Mean Curve (area is %.2f)' % back_mean_area)
    no_back_mean_fpr, no_back_mean_tpr = compute_mean_roc(no_back_fprs, no_back_tprs)
    no_back_mean_area = auc(no_back_mean_fpr, no_back_mean_tpr)
    plt.plot(no_back_mean_fpr, no_back_mean_tpr, lw=2, color='blue', linestyle='--', label='No Background Mean Curve (area is %.2f)' % no_back_mean_area)
    plt.xlim([-.05, 1.05])
    plt.ylim([-.05, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title(title)
    plt.legend(loc='lower right')
    plt.savefig(save_location)
    plt.close()

with open('all_assay_results/file_locations.pickle', 'rb') as f:
    locations = pickle.load(f)
    parser = argparse.ArgumentParser(description='Run method on dataset')
    parser.add_argument('prediction_method', help='Enter pwm for pwm, hmm for hmm, and quant for quantitative matrix')
    args = parser.parse_args()

    if args.prediction_method == 'pwm':
        method = PredictionMethod.PWM
    elif args.prediction_method == 'hmm':
        method = PredictionMethod.HMM
    else:
        print('Sorry, not a valid prediction method')
        sys.exit()

    assays = ['Kd_purified_direct_fluorescence', 'Kd_cellular_competitive_fluorescence', 'Kd_cellular_competitive_radioactivity', 'Kd_lysate_competitive _radioactivity', 'Kd_purified_competitive_fluorescence', 'Kd_purified_competitive_radioactivity', 'ic50_cellular_competitive_fluorescence', 'ic50_cellular_competitive_radioactivity', 'ic50_purified_competitive_fluorescence', 'ic50_purified_competitive_radioactivity']

    """For each assay, find the top 3 HLA's with the most data available. assay_hlas is of the form: {'assay':[(hla, filename, num_results, num_positive_results)...]}
    Make sure that there are at least 200 results, and at least 50 positive results.
"""  
    assay_hlas = dict()
    threshold = 500.0
    length = 9
    for assay in assays:
        assay_hlas[assay] = list()
        for hla, hla_assays in locations.items():
            if assay in hla_assays:
                filename = path.join('all_assay_results', hla_assays[assay])
                with open(filename, 'rb') as g:
                    assay_results = list(filter(lambda x: len(x[0]) == 9, pickle.load(g)))
                    num_results = len(assay_results)
                    num_positive_results = len(list(filter(lambda x: x[-1] <= threshold, assay_results)))
                    if num_results >= 200 and num_positive_results >= 50:
                        if len(assay_hlas[assay]) < 3:
                            assay_hlas[assay].append((hla, filename, num_results, num_positive_results))
                        elif num_results > min([x[2] for x in assay_hlas[assay]]):
                            #replace the HLA with the smallest number of results with this one.
                            value, index = min([(x[2], index) for index,x in enumerate(assay_hlas[assay])])
                            assay_hlas[assay][index] = (hla, filename, num_results, num_positive_results)

        for hla, filename, num_results, num_positive_results in assay_hlas[assay]:
            with open(filename, 'rb') as g:
                results = list(filter(lambda x: len(x[0]) == 9, pickle.load(g)))
                real_results = [(x[0], x[-1]) for x in results]
                title = 'Assay: {0}, HLA: {1}, \n threshold: {2}, Number of datapoints: {3}'.format(assay, hla, threshold, num_results)
                roc_curve_cross(real_results, threshold, title, 'background_comparison_' + args.prediction_method + '_' + assay + '_' + hla + '.png',  method, num_rows = 2)
                
                
