import requests
import json
import os
import pickle
import sys
from Bio.Alphabet import IUPAC
import argparse
import csv
import time

"""
The HLA allele and assay type are given as URLs such as: http://www.ontobee.org/ontology/MRO?iri=http://purl.obolibrary.org/obo/MRO_0001007

Specify one assay and HLA at a time.
"""
def get_id(hla_url, assay_url):
    print('hla url: {0}'.format(hla_url))
    print('assay url: {0}'.format(assay_url))
    params = {"structure_type":"linear_sequence","e_value":"exact","assay_t_cell_all":0,"assay_b_cell_all":0,"assay_mhc_all":1,"assay_mhc_ids[]":[assay_url],"mhc_class_type":"any_mhc","mhc_class_ids[]":[hla_url],"host_organism_type":"any_host","disease_type":"any_disease" ,"reference_type":"anyreference","count":"","page_num":"","sort_col":"","sort_dir":"","items_per_page":"","start_char":"","sort_col_type":"","list_type":"peptidic","search_type":"simple_search", "list_type":"mhcligand"}
    #We send in the params, and get back an id. We then query http://www.iedb.org/WorkRequestHandler.php with this ID. Do this until the result_code is not longer "PENDING"
    #If you don't do the json.dumps(params), this will not work!
    r = requests.post('http://www.iedb.org/WorkRequestHandler.php', data={'worker_type':2, 'params':json.dumps(params)})
    print('URL:')
    print(r.url)
    print('headers')
    print(r.request.headers)
    request_id = r.json()['id']
    return request_id


class IEDBRequestError(Exception):
    #Just pass in what was returned from the call
    def __init__(self, returned):
        self.message = 'There was an error with making the IEDB request. This was what was returned: {0}'.format(str(returned))

class CSVParseError(Exception):
    def __init__(self, error_message):
        self.message = error_message
        
class IEDBRequestStatus:
    def __init__(self, request_id):
        self._still_pending = True
        self._request_id = request_id
        self._error = False
    def setDownloadURL(self, url):
        self._still_pending = False
        self._download_url = url

    def isStillPending(self):
        return self._still_pending
    #Returns False if there is no download url
    def getDownloadURL(self):
        if self._still_pending:
            return False
        else:
            return self._download_url
    def getRequestID(self):
        return self._request_id

    def setErrorFlag(self):
        self._still_pending = False
        self._error = True

    def getErrorFlag(self):
        return self._error

    def setErrorInfo(self, info):
        self._error_info = info

    def getErrorInfo(self):
        return self._error_info
def check_request(iedb_request):
    r= requests.get('http://www.iedb.org/WorkRequestHandler.php', params={'id':iedb_request.getRequestID()})
    returned_data = r.json()
    if returned_data['result_code'] == 'SUCCESS':
        iedb_request.setDownloadURL('http://iedb.org' + returned_data['result_data'])
    elif returned_data['result_code'] != 'PENDING':
        iedb_request.setErrorFlag()
        iedb_request.setErrorInfo({'returned': returned_data, 'url':r.url})
    return iedb_request


"""
This function first checks which line (of the first 3) has 'Object Type', 'Description', 'Starting Position', 'Quantitative measurement', and 'Ending Position' in it. Then in takes the chunked up line, and prints it out to the user, asking for confirmation that the line it has chosen really does contain the field names. I'm doing this because the IEDB output seems to have two lines for field names, but that doesn't make any sense. 

Then, for each line past the fields line, we make sure that the Description field (which holds the peptide sequence) is a valid peptide, and that the Object Type is Linear Peptide, and that Ending Position - Starting Position + 1 = len(peptide). If any of these conditions fail, then we ask the user whether we should discard the data. 

Pass in the parsed lines of the CSV file. It will return a list of tuples of the form [(peptide, Kd)...]

"""
def parse_data(parsed_lines):
    print("I am about to extract the data from the CSV file. However, I will need you to confirm a few things with me along the way")
    i = 0
    fields_index = -1
    while i < 3:
        line = parsed_lines[i]
        print('line')
        print(line)
        if 'Object Type' in line and 'Description' in line and 'Starting Position' in line and 'Ending Position' in line and 'Quantitative measurement' in line:
            fields_index = i
            break
        i += 1
    if fields_index == -1:
        raise CSVParseError('We couldn\'t find the line with the field names in it. We searched the first 3 lines.')
    print('We selected these as the field names:')
    print(parsed_lines[fields_index])
#    keepGoing  = input('Is this correct? (yes/no) ')
 #   if keepGoing != 'yes':
 #       raise CSVParseError('User did not say that the field names were correct')

    data = list()
    #we add two, since line numbers are indexed at 1, and we are also at the line after the fields line.
    line_number = fields_index + 2
    for x in parsed_lines[(fields_index + 1)::]:
        entry = dict(zip(parsed_lines[fields_index], x))
        if 'Description' in entry and 'Starting Position' in entry and 'Ending Position' in entry and 'Object Type' in entry:
            peptide = entry['Description']
            object_type = entry['Object Type']
            assay_id = entry['MHC ligand ID']
            method = entry['Method/Technique']
            no_position = False
            try:
                starting_position = int(entry['Starting Position'])
            except ValueError:
                no_position = True
            try:
                ending_position = int(entry['Ending Position'])
            except ValueError:
                no_position = True
            position_criteria = no_position or ending_position - starting_position + 1 == len(peptide)
            try:
                kd = float(entry['Quantitative measurement'])
            except ValueError:
                print('could not convert quantitive measurement: {0}'.format(entry['Quantitative measurement']))
            else:
                if object_type == 'Linear peptide' and all([aa in IUPAC.protein.letters for aa in peptide]) and position_criteria:
                    print('added to data')
                    data.append((peptide, kd, assay_id, method))

        line_number += 1
    return data


"""

This returns a list of tuples of the from [(sequence, measurement)..]

Pass in the HLA and assay URL's
"""
def get_measurements(hla_url, assay_url):
    request_id  = get_id(hla_url, assay_url)
    request_status = IEDBRequestStatus(request_id)
    keepGoing = True
    while keepGoing:
        request_status = check_request(request_status)
        if request_status.getErrorFlag():
            raise IEDBRequestError(request_status.getErrorInfo())
        elif request_status.isStillPending() == False:
            keepGoing = False
            download_url = request_status.getDownloadURL()
            print('dowload url: {0}'.format(download_url))
        else:
            time.sleep(1)
    csv_request = requests.get(download_url)
    csv_text = csv_request.text
    lines = list(csv.reader(csv_text.split('\n')))
    return parse_data(lines)

parser = argparse.ArgumentParser(description='Get measurement data for assays and HLA alleles. This creates a directory for each HLA, and a file in each directory for each assay')
parser.add_argument('hla_pickle', help='Point us to a pickle that contains the HLA information, in the form: [(name, url)...], where name is the HLA name (such as HLA-A*01:01), and url is the OBI url')
parser.add_argument('assays_pickle', help='Point us to a pickle that contains the assay infomation, in the form: [(assay_name, url)...], where url is the OBI url')
args = parser.parse_args()
hla_pickle = args.hla_pickle
assays_pickle = args.assays_pickle
hla_list = False
if os.path.exists(hla_pickle):
    with open(hla_pickle, 'rb') as f:
        hla_list = pickle.load(f)
else:
    print('The hla_pickle file doesn\'t exist')
    sys.exit()

assay_list = False
if os.path.exists(assays_pickle):
    with open(assays_pickle, 'rb') as f:
        assay_list = pickle.load(f)
else:
    print('The assay_pickle file doesn\'t exist')
    sys.exit()

"""

"""    
file_locations = dict()
#When there's an exception with parsing the CSV stuff, add a dict with the fileds {'hla_name', 'hla_url', 'assay_name', 'assay_url', 'message'}
csv_exceptions = list()

for hla_name, hla_url in hla_list:
    file_locations[hla_name] = dict()
    os.mkdir(hla_name)
    for assay_name, assay_url in assay_list:
        try:
            measurements = get_measurements(hla_url, assay_url)
            file_location = os.path.join(hla_name, assay_name + '.pickle')
            with open(file_location, 'wb') as f:
                pickle.dump(measurements, f)
            file_locations[hla_name][assay_name] = file_location
        except CSVParseError as e:
            message = e.message
            print('THERE WAS A CSV PARSING ERROR')
            csv_exceptions.append({'hla_name': hla_name, 'hla_url':hla_url, 'assay_name':assay_name, 'assay_url':assay_url, 'message': message})


with open('file_locations.pickle', 'wb') as f:
    pickle.dump(file_locations, f)
