import pickle



with open('file_locations.pickle', 'rb') as f:
    locations = pickle.load(f)
    
#print('locations')
#print(locations)
for hla_name, assay_dict in locations.items():
    kd_location = assay_dict['Kd']
    with open(kd_location, 'rb') as g:
        kd = pickle.load(g)
    ic50_location = assay_dict['IC50']
    with open(ic50_location, 'rb') as g:
        ic50 = pickle.load(g)

    kd_sequence_set = set([sequence for sequence,measurement, assay_id, method, assay_group in kd])
    ic50_sequence_set = set([sequence for sequence, measurement, assay_id, method, assay_group in ic50])
    overlap = kd_sequence_set.intersection(ic50_sequence_set)
    if len(overlap) > 0:
        print('HLA: {0}'.format(hla_name))
        for sequence in overlap:
            print('Sequence: {0}'.format(sequence))
            for kd_sequence,measurement, assay_id, method, assay_group in kd:
                if kd_sequence == sequence:
                    print('Kd: {0}, {1}, {2}'.format(measurement, assay_id, method, assay_group))
            for ic_sequence, measurement, assay_id, method, assay_group in ic50:
                if ic_sequence == sequence:
                    print('IC50: {0}, {1}, {2}'.format(measurement, assay_id, method, assay_group))
