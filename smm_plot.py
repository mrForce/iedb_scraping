import matplotlib.pyplot as plt
import argparse
from sklearn.metrics import roc_curve, auc
import numpy as np
import pickle
from TrainingSystems import SMMUse
from os import path
from Bio.Alphabet import IUPAC
from itertools import cycle

from sklearn import svm, datasets

from sklearn.model_selection import KFold
from enum import Enum
import sys
class PredictionMethod(Enum):
    PWM = 1
    HMM = 2
    QuantMatrix = 3


"""
data should be of the form (sequence, measurement)

"""
def prediction_plot(data, title, save_location, sequence_length):
    #first, convert the measurements to True and False values based on the measurement threshold
    data = data
    sequences = np.array([sequence for sequence, measurement in data])
    print('sequences')
    print(sequences)
    measurements = np.array([measurement for sequence,measurement in data])
    print('measurements')
    print(measurements)
    cross_validation = KFold(n_splits=2)
    #    colors = cycle(['cyan', 'indigo', 'seagreen', 'yellow', 'blue'])
    colors = cycle(['cyan', 'green'])
    j = 1
    min_prediction = 10000
    max_prediction = -1000
    min_measurement = 10000
    max_measurement = -1000
    for (training_indices, test_indices), color in zip(cross_validation.split(sequences, measurements), colors):
        training_data_sequences = sequences[training_indices]
        training_data_measurements = measurements[training_indices]
        training_data = list(zip(training_data_sequences, training_data_measurements))
        print('about to train model {0}'.format(j))
        smm_model =  SMMUse(training_data, IUPAC.protein.letters, sequence_length)
        print('trained model {0}'.format(j))
        #predictions is a dictionary
        predictions = smm_model.predict([peptide for peptide in sequences[test_indices]])
        print('made predictions {0}'.format(j))
        #we turn it into an array. 
        test_predictions = np.array([predictions[peptide] for peptide in sequences[test_indices]])
        if min(test_predictions) < min_prediction:
            min_prediction = min(test_predictions)
        if max(test_predictions) > max_prediction:
            max_prediction = max(test_predictions)
        if min(measurements[test_indices]) < min_measurement:
            min_measurement = min(measurements[test_indices])
        if max(measurements[test_indices]) > max_measurement:
            max_measurement = max(measurements[test_indices])
        plt.scatter(test_predictions, measurements[test_indices], lw=2, color=color, label='Fold {0}'.format(j))
        j+= 1
    max_value = max(max_prediction, max_measurement)    
    plt.plot([0, max_value], [0, max_value], color='black', linestyle='--')
    plt.xlim([0, max_value])
    plt.gca().set_aspect('equal', adjustable='box')
    plt.ylim([0, max_value])
    plt.xlabel('Prediction')
    plt.ylabel('Measurement')
    plt.title(title)
    plt.legend(loc='lower right')
    plt.savefig(save_location)
    plt.close()

with open('all_assay_results/file_locations.pickle', 'rb') as f:
    locations = pickle.load(f)

    assays = ['Kd_purified_direct_fluorescence', 'Kd_cellular_competitive_fluorescence', 'Kd_cellular_competitive_radioactivity', 'Kd_lysate_competitive _radioactivity', 'Kd_purified_competitive_fluorescence', 'Kd_purified_competitive_radioactivity', 'ic50_cellular_competitive_fluorescence', 'ic50_cellular_competitive_radioactivity', 'ic50_purified_competitive_fluorescence', 'ic50_purified_competitive_radioactivity']

    """For each assay, find the top 3 HLA's with the most data available. assay_hlas is of the form: {'assay':[(hla, filename, num_results, num_positive_results)...]}
    Make sure that there are at least 200 results, and at least 50 positive results.
"""  
    assay_hlas = dict()
    threshold = 500
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
                title = 'Assay: {0}, HLA: {1}, \n Number of datapoints: {2}'.format(assay, hla, num_results)
                prediction_plot(real_results, title,  'smm_' + assay + '_' + hla + '.png', length)
