import matplotlib.pyplot as plt
import time
from TrainingSystems.smm import SMMTrainInput, TrainingDataType, SequenceDataType, DataPointType, SMMTrainOutput, SMMPredictor, SMMPredictInput, SMMPredictOutput, PredictType, parse, Greater, ThresholdType, Lesser
import tempfile
from scipy.stats import spearmanr
import subprocess
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

class SMMJob:
    """
    sequences is an np.array of strings.
    measurements is an np.array of measurements (floating point numbers) corresponding to the sequences in sequences
    train_indices is a list of indices, indicating which sequences/measurements are to be used for training.
    test_indices is a list of indices, indicating which sequences/measurements are to be used for testing
    alphabet is a string
    sequence_length is an integer
    predictor_save_location is a string that indicates where the predictor created by the SMM program should be saved.
    predictions_xml_save is a string that indicates where the predictions made by SMM -- that are formatted in XML -- should be saved.

    When we run a prediction, we should save of tuples of the form [(measurement, prediction),..]. 
    """
    def __init__(self, sequences, measurements, train_indices, test_indices, alphabet, sequence_length, predictor_save_location, predictions_xml_save, predictions_measurements_pickle_save):
        self.sequences = sequences
        self.measurements = measurements
        self.train_indices = train_indices
        self.test_indices = test_indices
        self.alphabet = alphabet
        self.sequence_length = sequence_length
        self.predictor_save_location = predictor_save_location
        self.predictions_xml_save = predictions_xml_save
        self.predictions_measurements_pickle_save = predictions_measurements_pickle_save
        self.training_running = 0
        self.predicting_running = 0 
        self.popen_process = None



    """
    This returns 0 if the job hasn't started to run.
    This returns 1 if the job is currently running.
    This returns 2 if the job finished
    """
    def get_training_running(self):
        return self.training_running
    
    def get_predicting_running(self):
        return self.predicting_running
    
    def set_popen_process(self, popen_process):
        self.popen_process = popen_process

    def get_popen_process(self):
        return self.popen_process
    

    def getTestIndices(self):
        return self.test_indices
    def getMeasurements(self):
        return self.measurements

    """
    This just sets the running flag to 2, and closes the temporary file we made.
    """
    def finish_training(self):
        assert(self.training_running == 1)
        self.popen_process.communicate()
        self.training_running = 2
        self.temp.close()
    def train(self):
        seq_data = SequenceDataType([DataPointType(Sequence = sequence, Measured=measurement) for sequence, measurement in zip(self.sequences[self.train_indices], self.measurements[self.train_indices])])
        training = TrainingDataType(Alphabet = ''.join(self.alphabet), SequenceLength = self.sequence_length, SequenceData = seq_data)
        smm_train = SMMTrainInput(OutputFile=self.predictor_save_location, TrainingData=training)
        self.temp = tempfile.NamedTemporaryFile('w+')
        with open(self.temp.name, 'w') as f:
            smm_train.export(f, 0)
        print('temp name in training: ' + self.temp.name)
        self.set_popen_process(subprocess.Popen(['smm', self.temp.name]))
        self.training_running = 1

    def predict(self):
        with open(self.predictor_save_location, 'r') as f:
            predictor = parse(f).get_SMMPredictor()
        predict_input = SMMPredictInput(OutputFile = self.predictions_xml_save, SMMPredictor = predictor, Predict=[PredictType(sequence) for sequence in self.sequences[self.test_indices]])
        self.temp = tempfile.NamedTemporaryFile()
        with open(self.temp.name, 'w') as f:
            predict_input.export(f, 0)
        print('temp name in predicting: ' + self.temp.name)
        self.set_popen_process(subprocess.Popen(['smm', self.temp.name]))
        self.predicting_running = 1
        
    def get_predictions(self):
        assert(self.predicting_running == 1)
        self.popen_process.communicate()
        self.predicting_running = 2
        with open(self.predictions_xml_save, 'r') as f:
            predict_output = parse(f).get_Predict()
            predictions = {pOutput.get_Sequence(): float(pOutput.get_Predictions()[0]) for pOutput in predict_output}
            print('predictions')
            print(predictions)
            print('save')
            print(self.predictions_xml_save)
            print('file input')
            with open(self.temp.name, 'r') as f:
                a = f.read()
                print(a)
            test_predictions = np.array([predictions[peptide] for peptide in self.sequences[self.test_indices]])
            """with open(self.predictions_measurements_pickle_save, 'w') as f:
                pickle.dump(np.column_stack((test_prediction, self.measurements[self.test_indices])), f)
            """
            return test_predictions
            
class SMMJobList:
    def __init__(self, title, plot_filename, zoomed_plot_filename, jobs):
        self.title = title
        self.plot_filename = plot_filename
        self.jobs = jobs
        self.zoomed_plot_filename = zoomed_plot_filename
    def getTitle(self):
        return self.title
    def getPlotFilename(self):
        return self.plot_filename
    def getZoomedPlotFilename(self):
        return self.zoomed_plot_filename
    def getJobs(self):
        return self.jobs
        
"""
data should be of the form (sequence, measurement)

"""
def prediction_plot(job_list, colors):
    plot_filename = job_list.getPlotFilename()
    zoomed_plot_filename = job_list.getZoomedPlotFilename()
    title = job_list.getTitle()
    
    min_prediction = 10000
    max_prediction = -1000
    min_measurement = 10000
    max_measurement = -1000
    j = 0
    for job, color in zip(job_list.getJobs(), colors):
        predictions = job.get_predictions()
        measurements = job.getMeasurements()[job.getTestIndices()]
        if min(predictions) < min_prediction:
            min_prediction = min(predictions)
        if max(predictions) > max_prediction:
            max_prediction = max(predictions)
        if min(measurements) < min_measurement:
            min_measurement = min(measurements)
        if max(measurements) > max_measurement:
            max_measurement = max(measurements)
        spearman, p_value = spearmanr(predictions, measurements)
        plt.scatter(predictions, measurements, lw=2, color=color, label='Fold {0}, (rho =  {1}, p =  {2})'.format(j, spearman, p_value))
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
    plt.savefig(plot_filename)
    plt.xlim([0, 5000])
    plt.ylim([0, 5000])
    plt.savefig(zoomed_plot_filename)
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
    sleep_time = 30
    max_num_jobs = 3
    colors = ['red', 'blue']
    job_lists = list()
    splits = 2
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
                real_results = [(x[0], x[-1]) for x in results][0:100]
                sequences = np.array([sequence for sequence, measurement in real_results])
                measurements = np.array([measurement for sequence, measurement in real_results])
                cross_validation = KFold(n_splits = splits)
                job_list = list()
                i = 0
                for training_indices, test_indices in cross_validation.split(sequences, measurements):
                    i += 1
                    name = 'smm_' + assay + '_' + hla
                    job = SMMJob(sequences, measurements, training_indices, test_indices, IUPAC.protein.letters, length, name + '_predictor_' + str(i) + '.xml', name + '_results_' + str(i) + '.xml', name + '_results_' + str(i)+ '.pickle')
                    job_list.append(job)
                    
                title = 'Assay: {0}, HLA: {1}, \n Number of datapoints: {2}'.format(assay, hla, num_results)
                filename = 'smm_' + assay + '_' + hla + '.png'
                zoomed_filename = 'smm_' + assay + '_' + hla + '_zoomed_.png'
                joblist_object = SMMJobList(title, filename, zoomed_filename, job_list)
                job_lists.append(joblist_object)
    
    
    """
    We want to store which jobs are running, and which ones are complete, in a matrix.
    Element at row r and column c, corresponding to job_lists[r].getJobs()[c], is:

    0 if job has not ran 
    1 if job is currently running (either training or predicting)
    2 if job has finished predicting
    3 if all jobs in that jobset are done predicting and making the plot.
    """
    jobs = np.zeros((len(job_lists), splits)).astype(int)
    i = 0
    while (0 in jobs) or (1 in jobs):
        """
        If there are less than max_num_jobs 1's in the matrix, and there are also some zeros, then we need to start some new jobs.
        """
        num_running_jobs = list(jobs.flatten()).count(1)
        if num_running_jobs < max_num_jobs and 0 in jobs:
            num_not_started = list(jobs.flatten()).count(0)
            if max_num_jobs - num_running_jobs > num_not_started:
                for row, column in np.transpose(np.nonzero(jobs == 0)):
                    job_lists[row].getJobs()[column].train()
                    jobs[row][column] = 1
            else:
                for row, column in np.transpose(np.nonzero(jobs == 0))[0:(max_num_jobs - num_running_jobs)]:
                    job_lists[row].getJobs()[column].train()
                    jobs[row][column] = 1
        """
        Now, go through all of the running jobs, and poll them to check if any have finished
        """
        for row, column in np.transpose(np.nonzero(jobs == 1)):
            job = job_lists[row].getJobs()[column]
            if job.get_popen_process().poll() != None:
                if job.get_training_running() == 1:
                    #If it was training, but now appears to be done.
                    job.finish_training()
                    job.predict()
                elif job.get_predicting_running() == 1:
                    #then done predicting
                    jobs[row][column] = 2
        """
        Now find any rows that are completely composed of 2s. That is, they're all done predicting.
        """
        for row_index in range(0, jobs.shape[0]):
            row = list(jobs[row_index][:])
            if row.count(2) == splits:
                prediction_plot(job_lists[row_index], colors)
                for column in range(0, splits):
                    #now set whole row to 3s, so we know it's done predicting
                    jobs[row_index][column] = 3
        #sleep for about 5 minutes
        print('jobs')
        print(jobs)
        time.sleep(sleep_time)
    
                    
