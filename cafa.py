import file_list as fl
import author_list as al
import taxon
import csv
import operator
import matplotlib.pyplot as plt
import numpy as np
import constant as const


class Cafa:

    def __init__(self):
        self.keyword_template = {
            'sequence alignment': 0, 'sequence-profile alignment': 0, 'profile-profile alignment': 0,
            'phylogeny': 0,'sequence properties': 0, 'physicochemical properties': 0, 'predicted properties': 0,
            'protein interactions': 0, 'gene expression': 0, 'mass spectrometry': 0, 'genetic interactions': 0,
            'protein structure': 0, 'literature': 0, 'genomic context': 0, 'synteny': 0, 'structure alignment': 0,
            'comparative model': 0, 'predicted protein structure': 0, 'de novo prediction': 0, 'machine learning': 0,
            'genome environment': 0, 'operon': 0, 'ortholog': 0, 'paralog': 0, 'homolog': 0, 'hidden Markov model': 0,
            'clinical data': 0, 'genetic data': 0, 'natural language processing': 0, 'other functional information': 0}

        self.data_root_directory = ''
        self.zip_data_directory = ''
        self.fmax_files_directory = ''
        self.csv_output_files_directory = ''
        self.taxon_name = 'None'
        self.set_path_cafa_team_files(const.DATA_ROOT_DIRECTORY)
        self.set_path_zipped_files(const.TAR_DATA_ROOT_DIRECTORY)
        self.set_path_fmax_files(const.FMAX_FILES_DIRECTORY)
        self.set_path_csv_out(const.CSV_OUTPUT_DIRECTORY)
        self.mode = 1   # Default mode.
        self.type = 1   # Default type.

        # initialize relative fmax scores dict and scores by equal weighting
        self.keyword_relative_fmax_score = {}
        self.keyword_equal_wts_score = {}
        self.keyword_relative_fmax_score = self.keyword_template.copy()
        self.keyword_equal_wts_score = self.keyword_template.copy()



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
        self.author_list = {}
        fl.create_file_list(author_file_list)
        al.create_author_list(self.author_list, author_file_list)

    # Pre-condition:
    # 1. Must have called create_author_kwd_taxon_dict() to create the author dictionary
    # 2. Must have called set_path_fmax_files
    def fmax_kwdscores_by_taxon(self, ontology, **kwargs):

        ontology = ontology.upper()

        # Set keywords scores dictionaries to zero for all keywords
        self.keyword_relative_fmax_score = self.keyword_template.copy()
        self.keyword_equal_wts_score = self.keyword_template.copy()
        self.fmax_sum = 0
        self.total_num_fmax_scores = 0

        # Default is type 1 mode 1.
        type = 1
        mode = 1

        if('type' in kwargs):
            type = kwargs['type']
        if('mode' in kwargs):
            mode = kwargs['mode']

        # Get either taxonName or taxonID from argument list
        if('taxonName' in kwargs and 'taxonID' in kwargs):
            print("Both taxonName and taxonID provided as arguments, only one is needed. Ignoring taxonID and using taxonName")
            self.taxon_name = kwargs['taxonName']
            self.taxon_id = taxon.taxon_ID_converter(self.taxon_name)
        elif('taxonName' in kwargs):
            self.taxon_name = kwargs['taxonName']
            self.taxon_id = taxon.taxon_ID_converter(self.taxon_name)
        elif('taxonID' in kwargs):
            self.taxon_name = taxon.taxon_name_converter(kwargs['taxonID'])
            self.taxon_id = kwargs['taxonID']
        else:
            raise Exception("'Either taxonName or taxonID needed as an argument.'")

        taxonID = taxon.taxon_ID_converter(self.taxon_name)


        #print("Ontology inside tabulator: "+o)
        file_name = self.fmax_files_directory + '/' + ontology.lower() + "_" + self.taxon_name + '_' + "type" + str(type) + '_' +\
                    'mode' + str(mode) + '_all_fmax_sheet.csv'


        with open(file_name, newline='') as csvfile:
            reader = csv.DictReader(csvfile)
            for i, row in enumerate(reader):
                if i > 1 and float(row['Coverage']) != 0.0:   # Skip if row is one of first 2 rows or coverage = 0
                    self.fmax_sum += float(row['F1-max'])
                    self.total_num_fmax_scores += 1

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
                                    self.keyword_relative_fmax_score[kwd] += float(row['F1-max'])/self.fmax_sum
                                    self.keyword_equal_wts_score[kwd] += 1.0/self.total_num_fmax_scores
                            else:
                                raise KeyError("Model " + model + " not found in author_list for taxonID: " +
                                               str(taxonID) + " Author: " + team)
                        else:
                            raise KeyError("TaxonID: " + str(taxonID) + "not found in author list for author: " + team)
                    else:
                        raise KeyError(team + " not found in author_list for taxonID: " + str(taxonID))
                        #print(str(i)+":"+team + " not found in author_list for taxonID: " + str(taxonID))

    # Pre-condition:
    # 1. Must have called create_author_kwd_taxon_dict() to create the author dictionary
    # 2. Must have called set_path_fmax_files
    # OLD CODE
    '''def fmax_kwdscores_by_taxonID(self, type, mode, taxonID):
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
            #print("Ontology inside tabulator: "+o)
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
        #for o in const.ONTOLOGY_LIST:
            #print("Ontology: "+o)
            #print("RELATIVE FMAX")
            #print(self.keyword_relative_fmax_score[o])
            #print("FMAX")
            #print(self.keyword_fmax_score[o])
            sorted_fmax_rel = sorted(self.keyword_relative_fmax_score[o].items(), key=operator.itemgetter(1), reverse=True)
            #sorted_fmax = sorted(self.keyword_fmax_score[o], key=operator.itemgetter(2), reverse=True)
            #print("SORTED FMAX")
            #print(sorted_fmax)
            #for item in sorted_fmax_rel:
                #print("%s ---> %.2f" % ( item[0], item[1]))

            #for item in sorted_fmax:
            #    print(item, self.keyword_fmax_score[o][item])
'''

    ''' Plots the relative frequency and the relative frequency weighted by fmax   
    '''
    def plot_ontology_fmax_results(self, ontology, **kwargs):

        type = 1
        mode = 1

        ontology = ontology.upper()

        if('type' in kwargs):
            type = kwargs['type']
        if('mode' in kwargs):
            mode = kwargs['mode']

#        if('taxonName' not in kwargs and 'taxonID' not in kwargs):
#            raise KeyError("Required argument missing. Need either taxonName or taxonID")
        # Get either taxonName or taxonID from argument list
        if('taxonName' in kwargs and 'taxonID' in kwargs):
            print("Both taxonName and taxonID provided as arguments, only one is needed. Ignoring taxonID and using taxonName")
            self.taxon_name = kwargs['taxonName']
            self.taxon_id = taxon.taxon_ID_converter(self.taxon_name)
        elif('taxonName' in kwargs):
            self.taxon_name = kwargs['taxonName']
            self.taxon_id = taxon.taxon_ID_converter(self.taxon_name)
        elif('taxonID' in kwargs):
            self.taxon_name = taxon.taxon_name_converter(kwargs['taxonID'])
            self.taxon_id = kwargs['taxonID']
        else:
            raise Exception("'Either taxonName or taxonID needed as an argument.'")

        plt.rc('xtick', labelsize=8)
        #fig, ax = plt.subplots(figsize=(10,2))
        #fig = plt.figure(figsize=(10,6))
        fig, ax = plt.subplots(figsize=(8, 6))
#        fig.suptitle("Ontology: "+ontology+"   Taxon: "+ self.taxon_name, fontsize=14)
        fig.canvas.set_window_title("Ontology: "+ontology+"   Taxon: "+ self.taxon_name)

        #ax = fig.add_axes([.1, .25, .8, .7])
        #ax.set_title(ontology + " Cumulative Relative Fmax Scores by Keyword for " + self.taxon_name + " Taxon")
        ax.set_ylabel("Relative Frequency")
        ax.set_xlabel("Keyword")

        bar_width = 0.4

        x_index = np.arange(0,len(self.keyword_relative_fmax_score))

        # Sorts the keywords of the Dictionary, self.keyword_relative_fmax_score, by the values of its keywords
        print("KEYWORD RELATIVE FMAX SCORE FOR ONTOLOGY "+ontology)
        print(self.keyword_relative_fmax_score)
        #sorted_fmax_rel = sorted(self.keyword_relative_fmax_score[ontology].items(), key=operator.itemgetter(1), reverse=True)
        sorted_scores = sorted(self.keyword_equal_wts_score.items(), key=operator.itemgetter(1), reverse=True)
        x_ticks = []#range(0,len(sorted_fmax_rel))
        y1 = []
        y2 = []
        for kwd, val in sorted_scores:
            x_ticks.append(kwd)
            y1.append(val)
            y2.append(self.keyword_relative_fmax_score[kwd])
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


    ''' Plots the relative frequency and the relative frequency weighted by fmax   
    '''
    '''def plot_ontology_fmax_results(self, ontology):

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
        print("KEYWORD RELATIVE FMAX SCORE FOR ONTOLOGY "+ontology)
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
        print(x_index+bar_width) '''


    ''' Precondition: Must have called self.fmax_kwdscores_by_taxonID() first.
    '''
    def csv_ontology_fmax_results(self, ontology, type, mode, taxonID):


        # Create CVS Header
        csvHeader = list(const.METHODOLOGY_KEYWORDS)
        csvHeader.insert(0, "Team")
        csvHeader.insert(1, "Taxon ID")
        csvHeader.insert(2, "Model Number")
        csvHeader.insert(3, "Fmax in "+ontology+"_"+"t"+str(type)+"_m"+str(mode))
        csvHeader.insert(4, "Equal weights")
        csvHeader.insert(5, "Relative Fmax")


        # Create csv input and output file names
        file_name = self.fmax_files_directory + '/' + ontology.lower() + "_" + self.taxon_name + '_' + "type" + str(type) + '_' +\
                    'mode' + str(mode) + '_all_fmax_sheet.csv'

        output_file_name = self.csv_output_files_directory + '/' + ontology.lower() + "_" + self.taxon_name + '_' + "type" + str(type) + '_' +\
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
                                out_dict["Fmax in " + ontology + "_" + "t" + str(type) + "_m" + str(mode)] = fmax
                                out_dict["Equal weights"] = 1.0 / self.total_num_fmax_scores
                                out_dict['Relative Fmax'] = float(fmax) / self.fmax_sum

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


    ''' Precondition: Must have called self.fmax_kwdscores_by_taxonID() first.
    '''
    '''def csv_ontology_fmax_results(self, type, mode, taxonID):

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

                        writer.writerow(out_dict)'''


    def create_taxonID_dictionary(self):
        #self.taxonID_counts = {}
        self.taxonID_counts = al.taxonID_counts(self.author_list).copy()

    def create_taxonID_list(self):
        self.taxonID_list = []

        for taxonID in self.taxonID_counts.keys():
            self.taxonID_list.append(taxonID)




    '''def fmax_kwdscores_by_ontology_taxon_list(self, type, mode, **kwargs):

        # Debug Code
        print("KWARGS")
        print(kwargs)

        print(kwargs['ontologies'])
        print(kwargs['bpo_taxons'])
        print(kwargs['cco_taxons'])
        print(kwargs['mfo_taxons'])

        taxa = {'10116': 'RAT', '9606': 'HUMAN', '3702': 'ARATH', '7955': 'DANRE', '44689': 'DICDI',
         '7227': 'DROME', '83333': 'ECOLI', '10090': 'MOUSE', '208963': 'PSEAE',
         '237561': 'CANAX', '559292': 'YEAST', '284812': 'SCHPO', '8355': 'XENLA', '224308': 'BACSU',
         '99287': 'SALTY', '243232': 'METJA', '321314': 'SALCH', '160488': 'PSEPK', '223283': 'PSESM',
         '85962': 'HELPY', '243273': 'MYCGE', '170187': 'STRPN', '273057': 'SULSO', 'all': 'all',
         'prokarya': 'prokarya', 'eukarya': 'eukarya'}

        # Debug Code
#        print("TAXA TABLE")
#        for k in taxa.keys():
#            print("taxa: " + k + "\tconverted: " + taxon.taxon_name_converter(k) + "\tunconverted back :" + taxon.taxon_ID_converter(taxon.taxon_name_converter(k)))

        #taxon.taxon_name_converter()
        # End Debug Code

        # Create taxon list for each ontology using arguments from **kwargs.
        taxon_list = {}
        for ontology in kwargs['ontologies']:
            taxon_list[ontology] = kwargs[ontology.lower()+"_taxons"]

        # Debug Code
        print("\n\nnewly created taxon list")
        for o in kwargs['ontologies']:
            print(o)
            print(taxon_list[o])
        # End Debug Code

        keyword_score_template = {
            'sequence alignment': 0, 'sequence-profile alignment': 0, 'profile-profile alignment': 0,
            'phylogeny': 0,'sequence properties': 0, 'physicochemical properties': 0, 'predicted properties': 0,
            'protein interactions': 0, 'gene expression': 0, 'mass spectrometry': 0, 'genetic interactions': 0,
            'protein structure': 0, 'literature': 0, 'genomic context': 0, 'synteny': 0, 'structure alignment': 0,
            'comparative model': 0, 'predicted protein structure': 0, 'de novo prediction': 0, 'machine learning': 0,
            'genome environment': 0, 'operon': 0, 'ortholog': 0, 'paralog': 0, 'homolog': 0, 'hidden Markov model': 0,
            'clinical data': 0, 'genetic data': 0, 'natural language processing': 0, 'other functional information': 0}

#        self.create_taxonID_dictionary()
#        self.create_taxonID_list()

        # Initialize the dictionaries keyword_relative_fmax_score and keyword_equal_wts_fmax_score
        keyword_relative_fmax_score = {}
        keyword_equal_wts_score = {}
#        for o in const.ONTOLOGY_LIST:
        for o in kwargs['ontologies']:
            keyword_relative_fmax_score[o] = keyword_score_template.copy()
            keyword_equal_wts_score[o] = keyword_score_template.copy()

        # Add entry to keyword_relative_fmax_score and keyword_equal_wts_fmax_score which will be for ALL
        # ontologies combined.
        keyword_relative_fmax_score['ALL'] = keyword_score_template.copy()
        keyword_equal_wts_score['ALL'] = keyword_score_template.copy()

        fmax_sum_all_ontologies = 0
        total_num_fmax_scores_all_ontologies = 0
        fmax_sum = {}
        total_num_fmax_scores = {}

        #for o in const.ONTOLOGY_LIST:
        for o in kwargs['ontologies']:

            fmax_sum[o] = 0
            total_num_fmax_scores[o] = 0

            for taxon_name in taxon_list[o]:

                #taxon_name = taxon.taxon_name_converter(taxonID)
                file_name = self.fmax_files_directory + '/' + o.lower() + "_" + taxon_name + '_' + "type" + str(type) + '_' +\
                            'mode' + str(mode) + '_all_fmax_sheet.csv'

                with open(file_name, newline='') as csvfile:
                    reader = csv.DictReader(csvfile)
                    for i, row in enumerate(reader):
                        if i > 1 and float(row['Coverage']) != 0.0:   # Skip if row is one of first 2 rows or coverage = 0
                            fmax_sum[o] += float(row['F1-max'])
                            total_num_fmax_scores[o] += 1
                            fmax_sum_all_ontologies += float(row['F1-max'])
                            total_num_fmax_scores_all_ontologies += 1

        print('RESULTS')
        for o in kwargs['ontologies']:
            print(fmax_sum[o])
            print(total_num_fmax_scores)
        print(fmax_sum_all_ontologies)
        print(total_num_fmax_scores_all_ontologies)

        #for o in const.ONTOLOGY_LIST:
        for o in kwargs['ontologies']:

            for taxon_name in taxon_list[o]:

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
                                            keyword_relative_fmax_score[o][kwd] += float(row['F1-max']) / fmax_sum[o]
                                            keyword_equal_wts_score[o][kwd] += 1.0 / total_num_fmax_scores[o]
                                            keyword_relative_fmax_score['ALL'][kwd] += float(row['F1-max']) / fmax_sum_all_ontologies
                                            keyword_equal_wts_score['ALL'][kwd] += 1.0 / total_num_fmax_scores_all_ontologies

                                    else:
                                        raise KeyError("Model " + model + " not found in author_list for taxonID: " +
                                                       str(taxonID) + " Author: " + team)
                                else:
                                    raise KeyError("TaxonID: " + str(taxonID) + "not found in author list for author: " + team)
                            else:
                                raise KeyError(team + " not found in author_list for taxonID: " + str(taxonID))
                                #print(str(i)+":"+team + " not found in author_list for taxonID: " + str(taxonID))


        # Debug Code
        #for o in const.ONTOLOGY_LIST:
            #print("Ontology: "+o)
            #print("RELATIVE FMAX")
            #print(self.keyword_relative_fmax_score[o])
            #print("FMAX")
            #print(self.keyword_fmax_score[o])
            sorted_fmax_rel = sorted(self.keyword_relative_fmax_score[o].items(), key=operator.itemgetter(1), reverse=True)
            #sorted_fmax = sorted(self.keyword_fmax_score[o], key=operator.itemgetter(2), reverse=True)
            #print("SORTED FMAX")
            #print(sorted_fmax)
            #for item in sorted_fmax_rel:
                #print("%s ---> %.2f" % ( item[0], item[1]))

            #for item in sorted_fmax:
            #    print(item, self.keyword_fmax_score[o][item]) '''


