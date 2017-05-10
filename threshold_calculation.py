import sys
import pickle
import os.path
import argparse

parser = argparse.ArgumentParser(description='This script takes in the results from MHC-ligand experiments, and t-cell activation experiments, and finds results between the two groups that used the same peptide and HLA allele. Suppose we have experiments that show peptide-HLA half life. Let set A be the set of half lives such that a t-cell activation occured for the same peptide-HLA combination. This script should calculate a threshold such that some percentage (specified by the user) of these half lives have a half life greater than or equal to this half life')
parser.add_argument('tcell_location_pickle', help='When you pulled t cell results from IEDB, you should have gotten a file_locations.pickle file that specifies where each of the tcell results for each HLA is located')
parser.add_argument('ligand_results_pickle', help='Pickle describes where the MHC-ligand assay results are stored')
parser.add_argument('assays_list_pickle', help='a pickle containing a list of tuples, of the form [(assay_name, url, ceiling)..] ceiling is True if the threshold should be a maximum (in the case of IC50 or Kd), False if the threshold should be a minimum.')
parser.add_argument('fractional_threshold', help='Put in 0.9 if you want 90% of the MHC-peptide measurements that are paired with positive t cell results to be contained by the threshold', type=float)
parser.add_argument('min_peptide_overlap', help='The minimum number of peptides that overlap between the MHC-ligand dataset and the positive t-cell dataset needed to be able to compute a threshold', type=int)
parser.add_argument('output_file', help='Where we output the thresholds')

args = parser.parse_args()
tcell_pickle_location = args.tcell_location_pickle
if not os.path.isfile(tcell_pickle_location):
    print('The tcell_location_pickle argument does not point to a valid file')
    sys.exit()

ligand_pickle = args.ligand_results_pickle
if not os.path.isfile(ligand_pickle):
    print('The ligand_result_pickle argument does not point to a valid file')
    sys.exit()

assays_pickle = args.assays_list_pickle
if not os.path.isfile(assays_pickle):
    print('The assays_list_pickle argument does not point to a valid file')
    sys.exit()

fractional_threshold = args.fractional_threshold
if fractional_threshold < 0.0:
    print('Fractional threshold must be greater than or equal to zero')
    sys.exit()
if fractional_threshold > 1.0:
    print('Fractional threshold must be less than or equal to one')
    sys.exit()

min_peptide_overlap = args.min_peptide_overlap
if min_peptide_overlap <= 0:
    print('min_peptide_overlap must be greater than zero')
    sys.exit()

assays = list()
with open(assays_pickle, 'rb') as f:
    assays = pickle.load(f)

tcell_locations = dict()
with open(tcell_pickle_location, 'rb') as f:
    tcell_locations = pickle.load(f)

tcell_location_head = os.path.split(tcell_pickle_location)[0]
ligand_results_locations = dict()
with open(ligand_pickle, 'rb') as f:
    ligand_results_locations = pickle.load(f)


ligand_pickle_head = os.path.split(ligand_pickle)[0]
#map assays to their threshold
thresholds = dict()
for assay_name, url, ceiling in assays:
    #we store the MHC-ligand results, for MHC-peptide pairs that match with a positive t cell event.
    measures = list()
    negative_measures = list()
    overlap_peptides = list()
    """ligand_results is a dict with the following form:
    {'HLA':{'assay_name':'location'...}...}
    """
    for hla, ligand_assays in ligand_results_locations.items():
        if hla in tcell_locations and assay_name in ligand_assays:
            if 'ic50' in assay_name:
                print('hla: {0}'.format(hla))
            with open(os.path.join(tcell_location_head, tcell_locations[hla]['tcells']), 'rb') as f:
                #of the form:  [['SADNNNSEY', 'Linear peptide', '142', 'in vitro assay', 'cytotoxicity', 'Positive-Low']...]
                tcell_results = pickle.load(f)
                #We want to get the peptides that, in combination with the given HLA allele, yielded a positive tcell result
                positive_peptides = set([x[0] for x in filter(lambda y: 'positive' in y[5].lower(), tcell_results)])
                #These peptides yielded ONLY negative results. 
                negative_peptides = set([x[0] for x in filter(lambda y: y[0] not in positive_peptides, tcell_results)])
                with open(os.path.join(ligand_pickle_head, ligand_assays[assay_name]), 'rb') as g:
                    ligand_results= pickle.load(g)
                    for ligand_result in ligand_results:
                        peptide = ligand_result[0]
                        measure = ligand_result[-1]
                        if peptide in positive_peptides:
                            if peptide not in overlap_peptides:
                                overlap_peptides.append(peptide)
                            measures.append(measure)
                        elif peptide in negative_peptides:
                            negative_measures.append(measure)
    if len(overlap_peptides) >= min_peptide_overlap:
        #if the threshold should be a ceiling, then sort in ascending order
        #if the threshold should be a floor, then sort in descending order
        measures.sort(reverse=(not ceiling))
        print('assay name')
        print(assay_name)
        print('measures')
        print(measures)
        print('negative measures')
        print(negative_measures)

        threshold = measures[int(fractional_threshold*len(measures)) - 1]
        print('threshold')
        print(threshold)
        #thresholds[assay_name] = threshold
        #count the number of measurements of pMHC binding (that yielded a negative when T-cell activation was measured) that would be classified as positive by the threshold


        num_false_positives = len(list(filter(lambda x: x <= threshold if ceiling else x >= threshold, negative_measures)))
        thresholds[assay_name] = (threshold, 1.0*num_false_positives/len(negative_measures), len(negative_measures))

with open(args.output_file, 'w') as f:
    f.write('fractional threshold: {0}\n'.format(fractional_threshold))
    for assay_name, (threshold, false_positive_rate, negative_measures) in thresholds.items():
        f.write('{0}: {1}, {2}, {3}\n'.format(assay_name, threshold, false_positive_rate, negative_measures))
