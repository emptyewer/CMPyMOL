import re
import copy
import math
import numpy as np
import functions.spinbar as spinbar

class pmap():

    index_ca = 0
    index_cb = 0

    def __init__(self, pdb_path, cutoff, use_ca, parent):
        self.aminoacids = ['GLU', 'ASP', 'LYS', 'ARG', 'HIS', 'GLN', 'PRO', 'ASN', 'ALA', 'THR', 'SER', 'VAL', 'GLY',
                           'MET', 'CYS', 'ILE', 'LEU', 'TYR', 'PHE', 'TRP']
        self.parent = parent
        self.pdb_path = pdb_path
        self.cutoff = cutoff
        self.use_ca = use_ca
        self.model_count = 0
        self.ca_coordinates = {}
        self.cb_coordinates = {}
        self.bfactor_cutoff = 25
        self.bfactors_ca = {}
        self.bfactors_cb = {}
        self.bfactors_present = False
        self.residue_names_ca = {}
        self.residue_numbers_ca = {}
        self.residue_names_cb = {}
        self.residue_numbers_cb = {}
        self.chain_names_ca = {}
        self.chain_names_cb = {}
        self.contact_maps_ca = {}
        self.contact_maps_cb = {}
        self.heat_maps = {}
        self.histogram_maps = {}
        self.variance_maps_ca = {}
        self.variance_maps_cb = {}
        self.parse_pdb()
        self.secondary_structure = []
        self.parent.append_status_text(">>> Calculating Maps...")
        _sbar = spinbar.SpinCursor(self.parent, msg="Calculating Maps... ")
        _sbar.start()
        self.calculate_contact_maps()
        self.calculate_variance_maps()
        _sbar.stop()



    def parse_pdb(self):
        global index_ca
        global index_cb

        model_match = re.compile(r'^MODEL')
        atom_match = re.compile(r'^ATOM')
        ca_match = re.compile(r'CA')
        cb_match = re.compile(r'CB')
        file_handle = open(self.pdb_path, 'r')


        def _initialize_dicts_for_model():
            global index_ca
            global index_cb
            index_ca = 0
            index_cb = 0
            self.model_count += 1
            self.bfactors_ca[self.model_count] = {}
            self.bfactors_cb[self.model_count] = {}
            self.ca_coordinates[self.model_count] = {}
            self.cb_coordinates[self.model_count] = {}
            self.residue_names_ca[self.model_count] = {}
            self.residue_numbers_ca[self.model_count] = {}
            self.residue_names_cb[self.model_count] = {}
            self.residue_numbers_cb[self.model_count] = {}
            self.chain_names_ca[self.model_count] = {}
            self.chain_names_cb[self.model_count] = {}
            self.variance_maps_ca[self.model_count] = {}
            self.variance_maps_cb[self.model_count] = {}
            self.histogram_maps[self.model_count] = {}

        for line in file_handle:
            if model_match.match(line):
                _initialize_dicts_for_model()

            if atom_match.match(line):
                if self.model_count == 0:
                    _initialize_dicts_for_model()

                line_split = line.split()

                shift = 0
                try:
                    resid = int(line_split[4])
                    shift = -1
                except ValueError:
                    resid = int(line_split[5])

                if ca_match.match(line_split[2]):
                    self.ca_coordinates[self.model_count][index_ca] = (float(line_split[6+shift]),
                                                                       float(line_split[7+shift]),
                                                                       float(line_split[8+shift]))
                    try:
                        self.bfactors_ca[self.model_count][index_ca] = float(line_split[10+shift])
                        if float(line_split[10 + shift]) > 0.0:
                            self.bfactors_present = True
                    except IndexError:
                        self.bfactors_ca[self.model_count][index_ca] = 0.0
                    self.residue_names_ca[self.model_count][index_ca] = line_split[3]
                    self.residue_numbers_ca[self.model_count][index_ca] = resid
                    if shift == 0:
                        self.chain_names_ca[self.model_count][index_ca] = line_split[4 + shift]
                    else:
                        self.chain_names_ca[self.model_count][index_ca] = "Z"
                    index_ca += 1
                    continue


                if cb_match.match(line_split[2]):
                    self.cb_coordinates[self.model_count][index_cb] = (float(line_split[6+shift]),
                                                                       float(line_split[7+shift]),
                                                                       float(line_split[8+shift]))
                    try:
                        self.bfactors_cb[self.model_count][index_cb] = float(line_split[10 + shift])
                        if float(line_split[10 + shift]) > 0.0:
                            self.bfactors_present = True
                    except IndexError:
                        self.bfactors_cb[self.model_count][index_cb] = 0.0

                    self.residue_names_cb[self.model_count][index_cb] = line_split[3]
                    self.residue_numbers_cb[self.model_count][index_cb] = resid
                    if shift == 0:
                        self.chain_names_cb[self.model_count][index_cb] = line_split[4 + shift]
                    else:
                        self.chain_names_cb[self.model_count][index_cb] = "Z"
                    index_cb += 1
                    continue

    def calculate_heatmap_histogram(self):
        if self.use_ca:
            contactmap = self.contact_maps_ca
            residue_names = self.residue_names_ca
            residue_numbers = self.residue_numbers_ca
        else:
            contactmap = self.contact_maps_cb
            residue_names = self.residue_names_cb
            residue_numbers = self.residue_numbers_cb

        for model in contactmap.keys():
            self.heat_maps[model] = {}
            self.histogram_maps[model] = np.zeros(len(residue_numbers[model]))
            self.heat_maps[model] = np.zeros((len(self.aminoacids), len(self.aminoacids)), np.int)
            cmap = copy.copy(contactmap[model])
            cmap[cmap > self.cutoff] = 0
            shape = cmap.shape

            for i in range(shape[0]):
                for j in range(i + 1, shape[1]):
                    if cmap[j][i] > 0:
                        self.histogram_maps[model][i] += 1
                        self.histogram_maps[model][j] += 1
                        index1 = self.aminoacids.index(residue_names[model][j])
                        index2 = self.aminoacids.index(residue_names[model][i])
                        self.heat_maps[model][index1][index2] += 1
                        if residue_names[model][j] != residue_names[model][i]:
                            self.heat_maps[model][index2][index1] += 1

    def calculate_variance_maps(self):
        def _calculate_variance_map(dictionary):
            mean = {}
            variance_map = {}
            mean[1] = np.zeros((dictionary[1].shape[0], dictionary[1].shape[1]))
            running_sum = np.zeros((dictionary[1].shape[0], dictionary[1].shape[1]))
            for model in range(1, self.model_count+1):
                running_sum += dictionary[model]
                mean[model] = running_sum/model

            sum_sub_sq = np.zeros((dictionary[1].shape[0], dictionary[1].shape[1]))
            for model in range(1, self.model_count+1):
                sub = dictionary[model] - mean[model]
                sum_sub_sq += np.square(sub)
                variance_map[model] = sum_sub_sq/model
            return variance_map

        self.variance_maps_ca = _calculate_variance_map(self.contact_maps_ca)
        self.variance_maps_cb = _calculate_variance_map(self.contact_maps_cb)


    def calculate_contact_maps(self):
        def _calc_dist(x, y):
            return math.sqrt((y[0] - x[0]) ** 2 + (y[1] - x[1]) ** 2 + (y[2] - x[2]) ** 2)

        def _calculate_contact_map(dictionary, coordinates):
            for model in range(1, self.model_count + 1):
                array_length = len(coordinates[model])
                keys = coordinates[model].keys()
                dist_array = np.zeros((array_length, array_length))
                for i in range(array_length):
                    for j in range(i + 1, array_length):
                        distance = _calc_dist(coordinates[model][keys[i]], coordinates[model][keys[j]])
                        # dist_array[i][j] = distance
                        dist_array[j][i] = distance
                dictionary[model] = dist_array

        _calculate_contact_map(self.contact_maps_ca, self.ca_coordinates)
        _calculate_contact_map(self.contact_maps_cb, self.cb_coordinates)






