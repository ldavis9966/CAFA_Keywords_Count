import file_list as fl
import author_list as al
import taxon
import csv
import operator
import matplotlib.pyplot as plt
import numpy as np
import constant as const

class cafa:

    #cafa
    #DATA_ROOT_DIRECTORY = 'data/raw_submission'
    #CSV_OUTPUT_DIRECTORY = 'csv_files'
    #TAR_DATA_ROOT_DIRECTORY = 'data'
    #TAR_FILE_NAME = 'raw_submission.tar.gz'

    def __init__(self):
        self.taxon_name = 'None'
        self.set_path_cafa_team_files(const.DATA_ROOT_DIRECTORY)
        self.set_path_zipped_files(const.TAR_DATA_ROOT_DIRECTORY)
        self.set_path_fmax_files(const.FMAX_FILES_DIRECTORY)
        self.set_path_csv_out(const.CSV_OUTPUT_DIRECTORY)


    def set_path_cafa_team_files(self, path):
        self.data_root_directory = path

    def set_path_zipped_files(self, path):
        self.zip_data_directory = path

    def set_path_fmax_files(self, path):
        self.fmax_files_directory = path

    def set_path_csv_out(self, path):
        self.csv_output_files_directory = path

    # Pre-condition:
    # Must have called self set_path_cafa_team_files('path') to ste path to cafa team files
    def create_author_kwd_taxon_dict(self):
        author_file_list = {}
        fl.create_file_list(author_file_list)
        self.author_list = {}
        al.create_author_list(self.author_list, author_file_list)

    # Pre-condition:
    # 1. Must have called create_author_kwd_taxon_dict() to create the author dictionary
    # 2. Must have called set_path_fmax_files
    def fmax_kwdscores_by_taxonID(self, type, mode, taxonID):
        self.keyword_fmax_score_template = {
            'sequence alignment': 0, 'sequence-profile alignment': 0, 'profile-profile alignment': 0,
            'phylogeny': 0,'sequence properties': 0, 'physicochemical properties': 0, 'predicted properties': 0,
            'protein interactions': 0, 'gene expression': 0, 'mass spectrometry': 0, 'genetic interactions': 0,
            'protein structure': 0, 'literature': 0, 'genomic context': 0, 'synteny': 0, 'structure alignment': 0,
            'comparative model': 0, 'predicted protein structure': 0, 'de novo prediction': 0, 'machine learning': 0,
            'genome environment': 0, 'operon': 0, 'ortholog': 0, 'paralog': 0, 'homolog': 0, 'hidden Markov model': 0,
            'clinical data': 0, 'genetic data': 0, 'natural language processing': 0, 'other functional information': 0}

        self.keyword_relative_fmax_score = {}
        self.keyword_fmax_score = {}

        for o in const.ONTOLOGY_LIST:
            self.keyword_relative_fmax_score[o] = self.keyword_fmax_score_template.copy()
            self.keyword_fmax_score[o] = self.keyword_fmax_score_template.copy()

        self.taxon_name = taxon.taxon_name_converter(taxonID)

        self.fmax_sum = {}
        self.total_num_fmax_scores = {}
        for o in const.ONTOLOGY_LIST:
            print("Ontology inside tabulator: "+o)
            file_name = self.fmax_files_directory + '/' + o.lower() + "_" + self.taxon_name + '_' + "type" + str(type) + '_' +\
                        'mode' + str(mode) + '_all_fmax_sheet.csv'

            self.fmax_sum[o] = 0
            self.total_num_fmax_scores[o] = 0
            with open(file_name, newline='') as csvfile:
                reader = csv.DictReader(csvfile)
                for i, row in enumerate(reader):
                    if i > 1 and float(row['Coverage']) != 0.0:   # Skip if row is one of first 2 rows or coverage = 0
                        self.fmax_sum[o] += float(row['F1-max'])
                        self.total_num_fmax_scores[o] += 1

            with open(file_name, newline='') as csvfile:
                reader = csv.DictReader(csvfile)
                for i, row in enumerate(reader):
                    if i > 1 and float(row['Coverage']) != 0.0:   # Skip if row is one of first 2 rows or coverage = 0
                        team = row['ID-model'][:-2].lower()
                        model = row['ID-model'][len(row['ID-model'])-1:]

                        if team in self.author_list:
                            if str(taxonID) in self.author_list[team]:
                                if model in self.author_list[team][str(taxonID)]:
                                    for kwd in self.author_list[team][str(taxonID)][model]:
                                        self.keyword_relative_fmax_score[o][kwd] += float(row['F1-max'])/self.fmax_sum[o]
                                        self.keyword_fmax_score[o][kwd] += 1.0/self.total_num_fmax_scores[o]
                                else:
                                    raise KeyError("Model " + model + " not found in author_list for taxonID: " +
                                                   str(taxonID) + " Author: " + team)
                            else:
                                raise KeyError("TaxonID: " + str(taxonID) + "not found in author list for author: " + team)
                        else:
                            raise KeyError(team + " not found in author_list for taxonID: " + str(taxonID))
                            #print(str(i)+":"+team + " not found in author_list for taxonID: " + str(taxonID))

        # Debug Code
        for o in const.ONTOLOGY_LIST:
            print("Ontology: "+o)
            print("RELATIVE FMAX")
            print(self.keyword_relative_fmax_score[o])
            print("FMAX")
            print(self.keyword_fmax_score[o])
            sorted_fmax_rel = sorted(self.keyword_relative_fmax_score[o].items(), key=operator.itemgetter(1), reverse=True)
            #sorted_fmax = sorted(self.keyword_fmax_score[o], key=operator.itemgetter(2), reverse=True)
            print("SORTED FMAX")
            #print(sorted_fmax)
            #for item in sorted_fmax_rel:
                #print("%s ---> %.2f" % ( item[0], item[1]))

            #for item in sorted_fmax:
            #    print(item, self.keyword_fmax_score[o][item])


    ''' Plots the relative frequency and the relative frequency weighted by fmax   
    '''
    def plot_ontology_fmax_results(self, ontology):

        ontology = ontology.upper()

        plt.rc('xtick', labelsize=8)
        #fig, ax = plt.subplots(figsize=(10,2))
        #fig = plt.figure(figsize=(10,6))
        fig, ax = plt.subplots(figsize=(8, 6))
        #ax = fig.add_axes([.1, .25, .8, .7])
        #ax.set_title(ontology + " Cumulative Relative Fmax Scores by Keyword for " + self.taxon_name + " Taxon")
        ax.set_ylabel("Relative Frequency")
        ax.set_xlabel("Keyword")

        bar_width = 0.4

        x_index = np.arange(0,len(self.keyword_relative_fmax_score[ontology]))

        # Sorts the keywords of the Dictionary, self.keyword_relative_fmax_score, by the values of its keywords
        print(self.keyword_relative_fmax_score[ontology])
        #sorted_fmax_rel = sorted(self.keyword_relative_fmax_score[ontology].items(), key=operator.itemgetter(1), reverse=True)
        sorted_fmax = sorted(self.keyword_fmax_score[ontology].items(), key=operator.itemgetter(1), reverse=True)
        x_ticks = []#range(0,len(sorted_fmax_rel))
        y1 = []
        y2 = []
        for kwd, val in sorted_fmax:
            x_ticks.append(kwd)
            y1.append(val)
            y2.append(self.keyword_relative_fmax_score[ontology][kwd])
        #plt.bar(x_index, y1, bar_width, color='#5C89C4', label="Relative FMax")
        #plt.bar(x_index+bar_width, y2, bar_width, color='#5DC687', label="FMax (Equally Weighted)")
        plt.bar(x_index, y1, bar_width, color='#5C89C4', label="Equal Weights")
        plt.bar(x_index+bar_width, y2, bar_width, color='#5DC687', label="Weighted by Fmax")

        index = np.arange(4)
        print(index)

        plt.legend()
        plt.xticks(x_index+bar_width/2, x_ticks)
        plt.xticks(rotation=90)
        plt.tight_layout()
        plt.show()
        print(x_index)
        print(x_index+bar_width)


    ''' Precondition: Must have called self.fmax_kwdscores_by_taxonID() first.
    '''
    def csv_ontology_fmax_results(self, type, mode, taxonID):

        for o in const.ONTOLOGY_LIST:
            # Create CVS Header
            csvHeader = list(const.METHODOLOGY_KEYWORDS)
            csvHeader.insert(0, "Team")
            csvHeader.insert(1, "Taxon ID")
            csvHeader.insert(2, "Model Number")
            csvHeader.insert(3, "Fmax in "+o+"_"+"t"+str(type)+"_m"+str(mode))
            csvHeader.insert(4, "Equal weights")
            csvHeader.insert(5, "Relative Fmax")


            # Create csv input and output file names
            file_name = self.fmax_files_directory + '/' + o.lower() + "_" + self.taxon_name + '_' + "type" + str(type) + '_' +\
                        'mode' + str(mode) + '_all_fmax_sheet.csv'

            output_file_name = self.csv_output_files_directory + '/' + o.lower() + "_" + self.taxon_name + '_' + "type" + str(type) + '_' +\
                        'mode' + str(mode) + '_all_fmax_sheet_with_weights.csv'

            # open both input and output csv file
            with open(file_name, newline='') as csvfile_in, open(output_file_name, 'w', newline='') as csvfile_out:

                reader = csv.DictReader(csvfile_in)
                writer = csv.DictWriter(csvfile_out, fieldnames=csvHeader)

                writer.writeheader()

                # Loop through each row in the input file and calculate the weights (equal and by fmax) then out
                # it to the new file.
                out_dict = {}
                for i, row in enumerate(reader):
                    if i > 1 and float(row['Coverage']) != 0.0:   # Skip if row is one of first 2 rows or coverage = 0

                        team = row['ID-model'][:-2].lower()
                        model = row['ID-model'][len(row['ID-model'])-1:]
                        fmax = row['F1-max']

                        out_dict.clear()
                        if team in self.author_list:
                            if str(taxonID) in self.author_list[team]:
                                if model in self.author_list[team][str(taxonID)]:

                                    out_dict.clear()
                                    out_dict['Team'] = team
                                    out_dict['Taxon ID'] = taxonID
                                    out_dict['Model Number'] = model
                                    out_dict["Fmax in " + o + "_" + "t" + str(type) + "_m" + str(mode)] = fmax
                                    out_dict["Equal weights"] = 1.0 / self.total_num_fmax_scores[o]
                                    out_dict['Relative Fmax'] = float(fmax) / self.fmax_sum[o]

                                    for kwd in self.author_list[team][str(taxonID)][model]:
                                        out_dict[kwd] = 1
                                else:
                                    raise KeyError("Model " + model + " not found in author_list for taxonID: " +
                                                   str(taxonID) + " Author: " + team)
                            else:
                                raise KeyError("TaxonID: " + str(taxonID) + "not found in author list for author: " + team)
                        else:
                            raise KeyError(team + " not found in author_list for taxonID: " + str(taxonID))

                        writer.writerow(out_dict)