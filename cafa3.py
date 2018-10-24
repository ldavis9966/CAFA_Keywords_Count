# import file_list as fl
# import author_list as al
import constant as const
import taxon

import csv
import operator
import matplotlib.pyplot as plt
import numpy as np
import re
import os


class Cafa3:

    def __init__(self):
        self.keyword_template = {
            'sequence alignment': 0, 'sequence-profile alignment': 0, 'profile-profile alignment': 0,
            'phylogeny': 0, 'sequence properties': 0, 'physicochemical properties': 0, 'predicted properties': 0,
            'protein interactions': 0, 'gene expression': 0, 'mass spectrometry': 0, 'genetic interactions': 0,
            'protein structure': 0, 'literature': 0, 'genomic context': 0, 'synteny': 0, 'structure alignment': 0,
            'comparative model': 0, 'predicted protein structure': 0, 'de novo prediction': 0, 'machine learning': 0,
            'genome environment': 0, 'operon': 0, 'ortholog': 0, 'paralog': 0, 'homolog': 0, 'hidden Markov model': 0,
            'clinical data': 0, 'genetic data': 0, 'natural language processing': 0, 'other functional information': 0}

        # File path vars for input and output.
        self.data_root_directory = ''
        self.zip_data_directory = ''
        self.fmax_files_directory = ''
        self.csv_output_files_directory = ''

        # Dictionary that contains all authors along with their taxon, model, and keywrod data. This drives everything
        # in the class
        self.author_list = {}

        # Taxon to be used for this instance. Note both name and ID must refer to the same taxon.
        self.taxon_name = 'None'  # Taxon name
        self.taxon_id = 0  # Taxon ID
        self.ontology = 'MFO'  # Default ontology

        # Taxa information derived from Cafa3 raw submission file and stored in self.author_list
        self.taxon_id_list = []  # List of all taxa in self.author_list
        self.taxon_id_counts = {}  # Number of occurrences for each taxon in self.author_list

        self.set_path_cafa_team_files(const.DATA_ROOT_DIRECTORY)
        self.set_path_zipped_files(const.TAR_DATA_ROOT_DIRECTORY)
        self.set_path_fmax_files(const.FMAX_FILES_DIRECTORY)
        self.set_path_csv_out(const.CSV_OUTPUT_DIRECTORY)
        self.bnchmrk_mode = 1   # Default mode.
        self.bnchmrk_type = 1   # Default type.

        # initialize relative fmax scores dict and scores by equal weighting
        self.keyword_relative_fmax_score = self.keyword_template.copy()
        self.keyword_equal_wts_score = self.keyword_template.copy()

        # Initialize raw counts dictionaries. One for total keywords counts and one for a particular
        # model (either model 1, 2, or 3)
        self.total_keyword_count = self.keyword_template.copy()
        self.model_keyword_count = self.keyword_template.copy()
        self.model = 0  # Set once a call to count_kwds() is made, used by subsequent call to plot_raw_keyword_counts()

        self.fmax_sum = 0
        self.total_num_fmax_scores = 0

        # Create author dictionary with taxon, model, and keyword data. Populates self.author_list.
        self.create_author_kwd_taxon_dict()

    def set_path_cafa_team_files(self, path):
        self.data_root_directory = path

    def set_path_zipped_files(self, path):
        self.zip_data_directory = path

    def set_path_fmax_files(self, path):
        self.fmax_files_directory = path

    def set_path_csv_out(self, path):
        self.csv_output_files_directory = path

    @staticmethod
    def create_author_list(authors_list, file_list):

        files_processed = 0

        for file in file_list:
            # Following code was originally designed to read in the cafa3 raw submission files using a tar file. The
            # performance was extremely slow (taking several hours to read all files) so not really useful
            # unless there is a way to boost performance.
            # if const.USE_ZIPPED_FILES:
            #    tarf = tarfile.open(const.TAR_DATA_ROOT_DIRECTORY + "/" + const.TAR_FILE_NAME)
            #    fh = tarf.extractfile(file_list[file]['path']+"/"+file)
            # else:
            #    fh = open(file_list[file]['path']+"/"+file, "r")
            fh = open(file_list[file]['path']+"/"+file, "r")

            # See if AUTHOR is in the first line. Extract the team name if so.
            cur_line = fh.readline()
            # if type(cur_line) is bytes:
            #    cur_line = cur_line.decode('utf-8')

            if 'AUTHOR' not in cur_line:
                print("\nWarning: AUTHOR not found in file.")
                Cafa3.printfile(file_list[file]['path'], file)
                break

            author_in_file = cur_line[len('AUTHOR '):].lower()
            author_in_file = author_in_file.strip()
            author_in_file = author_in_file.replace("\n", "")

            # See if MODEL is in the second line. Extract Model # if so.
            cur_line = fh.readline()
            # if type(cur_line) is bytes:
            #    cur_line = cur_line.decode('utf-8')

            if 'MODEL' not in cur_line:
                print("\nWarning: MODEL not found in file.")
                Cafa3.printfile(file_list[file]['path'], file)
                break

            model_in_file = cur_line[len('MODEL '):]
            model_in_file = model_in_file.replace("\n", "")
            model_in_file = model_in_file.strip(" ")

            file_str_split = file.split('_')
            author = file_str_split[0].strip().lower()
            model = file_str_split[1].strip()
            taxon_id = file_str_split[2].split('.')[0].strip(" ")

            # Check if model # and Author (Team Name) are the 1st two items in the file name per CAFA specifications
            if author_in_file.lower() != file_str_split[0].lower() or model_in_file != file_str_split[1]:
                print("\nWarning: filename conflicts with file data")
                print("\tFile: "+file_list[file]['path'] + "/" + file)
                # print("\tFilename string: ", file)
                # print("\tFile string parse:", file_str_split)
                if author_in_file.lower() != file_str_split[0].lower():
                    print("\t\tAuthor in filename does not match author shown in file")
                    print("\t\t\tAuthor (file name):", file_str_split[0])
                    print("\t\t\tAuthor (in file):", author_in_file)
                    # Assume the file name author is the correct author name
                if model_in_file != file_str_split[1]:
                    print("\t\tModel # in filename does not match model # shown in file")
                    print("\t\t\tModel (file name):", file_str_split[1])
                    print("\t\t\tModel (in file):", model_in_file)
                    # Assume the file name model # is the correct model number

            cur_line = fh.readline()
            # if type(cur_line) is bytes:
            #    cur_line = cur_line.decode('utf-8')

            if 'KEYWORDS' not in cur_line:
                print("\nWarning: KEYWORDS tag not found in file.")
                Cafa3.printfile(file_list[file]['path'], file)
                break

            # 3 part process: str.split() last part of current line(everything after the KEYWORDS tag) into a list using
            # comma delimiter to extract the keywords. Next remove all characters except upper & lower case alpha chars,
            # spaces, & hyphens from each keyword in the list. Finally str.strip() any leading or trailing spaces
            # from each keyword.
            keywords = cur_line[len('KEYWORDS '):]
            keywords = keywords.split(",")
            regex = re.compile('[^a-zA-Z -]')
            for i, kwd in enumerate(keywords):
                keywords[i] = regex.sub('', kwd)
                keywords[i] = keywords[i].strip(" ")

            # Check if the keywords in this file are one of the Cafa3 accepted keywords
            for kwrd in keywords:
                if kwrd not in const.METHODOLOGY_KEYWORDS:
                    print("\nWarning: keyword", "\""+kwrd+"\"", "not an acceptable methodology keyword in file:")
                    Cafa3.printfile(file_list[file]['path'], file)
                    break

            # This code block checks if current author is not in the author list. If not present then create
            # a new nested dictionary of both the taxon_id and then another nested dictionary using the current model #.
            if author not in authors_list:
                authors_list[author] = {taxon_id: {model: {}}}
            # Otherwise the author already has an entry. Now check if the current taxon_id is present.
            # If not then create a double nested dictionary for the author using the current taxon_id and model #.
            elif taxon_id not in authors_list[author]:
                authors_list[author][taxon_id] = {model: {}}
            # Otherwise the author already has a taxon_id entry. Now check if the current model has already been
            # added for the current taxon_id. If not, then create a new dictionary for this taxon_id using
            # the current model #.
            elif model not in authors_list[author][taxon_id]:
                authors_list[author][taxon_id][model] = {}
            # If all the prior checks fail, then there must be a double entry of the model # for current taxon_id.
            # So report error. In this case the current Model # will overwrite the prior model # for this taxon_id.
            # Probably should handle this better in the future.
            else:
                print("\nWarning: duplicate Model # found for Taxon ID")
                print("\tTaxonID:", taxon_id)
                print("\tModel #:", model)
                print("\tConflicting file:")
                print("\t\t", file)

            # For every keyword in this current text file, increment the keyword count and store it in the
            # authors_list for this particular author. Also update the total global keyword count.
            for kwrd in keywords:
                authors_list[author][taxon_id][model][kwrd] = 1

            # Close current file and increase counter
            fh.close()
            files_processed += 1
        return files_processed

    @staticmethod
    def printfile(path, filename):
        print("File Path + Name:", path+"/"+filename)
        return

    # Returns a dictionary where the number of occurences for each taxon_id is of the form
    # the {Key: value} = {taxon_id: count}
    @staticmethod
    def get_taxon_id_counts(author_list):
        taxon_id_counts = {}
        for author in author_list:
            for taxon_id in author_list[author]:
                if taxon_id not in taxon_id_counts:
                    taxon_id_counts[taxon_id] = 1
                else:
                    taxon_id_counts[taxon_id] += 1
        return taxon_id_counts

    # Directory walk the data folders (from either a zipped file or unzipped file) and place all the CAFA filenames in a
    # Python nested dictionary with structure as follows { filename : {'path': filepath} }.
    @staticmethod
    def create_file_list(file_list):
        file_count = 0

        # Following code was originally designed to read in the cafa3 raw submission files using a tar file. The
        # performance was extremely slow (taking several hours just to read in all the file names) so not really useful
        # unless there is a way to boost performance.
        # if const.USE_ZIPPED_FILES:  # use zipped files
        #    tarf = tarfile.open(constant.TAR_DATA_ROOT_DIRECTORY + "/" + constant.TAR_FILE_NAME)
        #        tarf = tarfile.open("data.tar.gz")
        #    for tarinfo in tarf:
        #        if tarinfo.isreg():                 # Check tarinfo is a file, and not a directory
        #            tmp = tarinfo.name.split("/")   # Split out the path+filename strings into a tuple using '/'
        #            filename = tmp[len(tmp) - 1]    # Get filename from last tuple
        #            del tmp[-1]                     # last entry from tuple which is the filename.
        #            path = ""                       # now rebuild the path from the tuple.
        #            for directory in tmp:
        #                path = path + directory + "/"
        #            path = path.strip("/")
        #            file_list[filename] = {'path': path}
        #            file_count += 1
        #            sb.print(2740, 20, file_count)
        #    tarf.close()
        # else:  # Use unzipped files
        #    for dirpath, dirnames, filenames in os.walk(constant.DATA_ROOT_DIRECTORY):
        #        for file in filenames:
        #            file_list[file] = {'path': dirpath}
        #            file_count += 1

        for dirpath, dirnames, filenames in os.walk(const.DATA_ROOT_DIRECTORY):
            for file in filenames:
                file_list[file] = {'path': dirpath}
                file_count += 1

        return file_count

    # Pre-condition:
    # Must have called self set_path_cafa_team_files('path') to ste path to cafa team files
    def create_author_kwd_taxon_dict(self):
        author_file_list = {}
        self.author_list = {}
        Cafa3.create_file_list(author_file_list)
        Cafa3.create_author_list(self.author_list, author_file_list)

    # Required argument
    # ontology - Must be either BPO, MFO, or CCO
    #
    # Optional arguments:
    # type - Either type 1 or 2. If not set the class default is 1
    # mode - Either mode 1 or 2.
    #
    #
    # Pre-condition:
    # 1. Must have called create_author_kwd_taxon_dict() to create the author dictionary
    # 2. Must have called set_path_fmax_files
    def kwdscores_by_taxon(self, ontology, **kwargs):

        ontology = ontology.upper()
        self.ontology = ontology

        # Set keywords scores dictionaries to zero for all keywords
        self.keyword_relative_fmax_score = self.keyword_template.copy()
        self.keyword_equal_wts_score = self.keyword_template.copy()
        self.fmax_sum = 0
        self.total_num_fmax_scores = 0

        # Default is type 1 mode 1.
        self.bnchmrk_type = 1
        self.bnchmrk_mode = 1

        if 'type' in kwargs:
            if kwargs['type'] != 1 and kwargs['type'] != 2:
                raise Exception("Invalid benchmark type given. Expect either 1 or 2.")
            self.bnchmrk_type = kwargs['type']
        if'mode' in kwargs:
            if kwargs['mode'] != 1 and kwargs['type'] != 2:
                raise Exception("Invalid benchmark mode given. Expect either 1 or 2.")
            self.bnchmrk_mode = kwargs['mode']

        # Get either taxonName or taxon_id from argument list
        if'taxonName' in kwargs and 'taxonID' in kwargs:
            print("Both taxonName and taxonID provided as arguments, only one is needed. Ignoring taxonID and using" +
                  " taxonName")
            self.taxon_name = kwargs['taxonName']
            self.taxon_id = taxon.taxon_name_to_id(self.taxon_name)
        elif 'taxonName' in kwargs:
            self.taxon_name = kwargs['taxonName']
            self.taxon_id = taxon.taxon_name_to_id(self.taxon_name)
        elif 'taxonID' in kwargs:
            self.taxon_name = taxon.taxon_id_to_name(kwargs['taxonID'])
            self.taxon_id = kwargs['taxonID']
        else:
            raise Exception("'Either taxonName or taxonID needed as an argument.'")

        # taxon_id = taxon.taxon_ID_converter(self.taxon_name)

        # print("Ontology inside tabulator: "+o)
#        if ontology.lower == 'all':
#            file_name = self.fmax_files_directory + '/' + ontology.lower() + "_all_type" + str(bnchmrk_type) + '_' +\
#                    'mode' + str(bnchmrk_mode) + '_all_fmax_sheet.csv'
#        else:
#        file_name = self.fmax_files_directory + '/' + ontology.lower() + "_" + self.taxon_name + '_' + "type" +\
#            str(self.bnchmrk_type) + '_' + 'mode' + str(self.bnchmrk_mode) + '_all_fmax_sheet.csv'

        # For Linux path
        file_name = self.fmax_files_directory + '/' + ontology.lower() + "_" + self.taxon_name + '_' + "type" +\
            str(self.bnchmrk_type) + '_' + 'mode' + str(self.bnchmrk_mode) + '_all_fmax_sheet.csv'

        # For Windows 10 path
        # file_name = self.fmax_files_directory + '\\' + ontology.lower() + "_" + self.taxon_name + '_' + "type" +\
        #    str(self.bnchmrk_type) + '_' + 'mode' + str(self.bnchmrk_mode) + '_all_fmax_sheet.csv'


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

                    if self.taxon_name == 'all':
                        # print("A \'"+self.taxon_name+"\'")
                        if team in self.author_list:
                            for taxon_id in self.author_list[team]:
                                # print("Author: "+team+"\ttaxon_id: "+taxon_id)
                                # for i, key in enumerate(self.author_list[team][str(taxon_id)]):
                                    # print("Key "+str(i)+" : "+key)

                                if model in self.author_list[team][str(taxon_id)]:
                                    for kwd in self.author_list[team][str(taxon_id)][model]:
                                        self.keyword_relative_fmax_score[kwd] += float(row['F1-max'])/self.fmax_sum
                                        self.keyword_equal_wts_score[kwd] += 1.0/self.total_num_fmax_scores
                        else:
                            raise KeyError(team + " not found in author_list for taxonID: " + str(taxon_id))

                    else:
                        # print("B \'"+self.taxonName+"\'")
                        if team in self.author_list:
                            if str(self.taxon_id) in self.author_list[team]:
                                if model in self.author_list[team][str(self.taxon_id)]:
                                    for kwd in self.author_list[team][str(self.taxon_id)][model]:
                                        self.keyword_relative_fmax_score[kwd] += float(row['F1-max'])/self.fmax_sum
                                        self.keyword_equal_wts_score[kwd] += 1.0/self.total_num_fmax_scores
                                else:
                                    raise KeyError("Model " + model + " not found in author_list for taxonID: " +
                                                   str(self.taxon_id) + " Author: " + team)
                            else:
                                raise KeyError("TaxonID: " + str(self.taxon_id) + "not found in author list for author: "
                                               + team)
                        else:
                            raise KeyError(team + " not found in author_list for taxonID: " + str(self.taxon_id))
                        # print(str(i)+":"+team + " not found in author_list for taxonID: " + str(taxon_id))

    ''' Plots the relative frequency and the raw counts weighted by fmax   
    '''
#    def plot_ontology_score_results(self, ontology, **kwargs):
    def plot_ontology_score_results(self, width, height, x_font_size, y_font_size, label_font_size, bar1_hex_color, bar2_hex_color):

#        self.bnchmrk_type = 1
#        self.bnchmrk_mode = 1

        ontology = self.ontology.upper()

#        if 'type' in kwargs:
#            self.bnchmrk_type = kwargs['type']
#        if'mode' in kwargs:
#            self.bnchmrk_mode = kwargs['mode']

        print("Mode: "+str(self.bnchmrk_mode))
        print("Type: "+str(self.bnchmrk_type))

#        if('taxonName' not in kwargs and 'taxonID' not in kwargs):
#            raise KeyError("Required argument missing. Need either taxonName or taxonID")
        # Get either taxonName or taxonID from argument list
#        if 'taxonName' in kwargs and 'taxonID' in kwargs:
#            print("Both taxonName and taxonID provided as arguments, only one is needed. Ignoring taxonID and " +
#                  "using taxonName")
#            self.taxon_name = kwargs['taxonName']
#            self.taxon_id = taxon.taxon_name_to_id(self.taxon_name)
#        elif 'taxonName' in kwargs:
#            self.taxon_name = kwargs['taxonName']
#            self.taxon_id = taxon.taxon_name_to_id(self.taxon_name)
#        elif 'taxonID' in kwargs:
#            self.taxon_name = taxon.taxon_id_to_name(kwargs['taxonID'])
#            self.taxon_id = kwargs['taxonID']
#        else:
#            raise Exception("'Either taxonName or taxonID needed as an argument.'")

        plt.rc('xtick', labelsize=x_font_size)
        plt.rc('ytick', labelsize=y_font_size)
        plt.rcParams.update({'font.size': label_font_size})

        plt.rc('xtick', labelsize=8)
        # fig, ax = plt.subplots(figsize=(10,2))
        # fig = plt.figure(figsize=(10,6))
        # fig, ax = plt.subplots(figsize=(8, 6))
        fig, ax = plt.subplots(figsize=(width, height))
#        fig.suptitle("Ontology: "+ontology+"   Taxon: "+ self.taxon_name, fontsize=14)
        fig.canvas.set_window_title("Ontology: " + ontology + "\tTaxon: " + self.taxon_name+"\tType: " +
                                    str(self.bnchmrk_type) + "\tMode: " + str(self.bnchmrk_mode))

        # ax = fig.add_axes([.1, .25, .8, .7])
        # ax.set_title(ontology + " Cumulative Relative Fmax Scores by Keyword for " + self.taxon_name + " Taxon")
        ax.set_ylabel("Relative Frequency")
        ax.set_xlabel("Keyword")

        bar_width = 0.4

        x_index = np.arange(0, len(self.keyword_relative_fmax_score))

        # Sorts the keywords of the Dictionary, self.keyword_relative_fmax_score, by the values of its keywords
        print("KEYWORD RELATIVE FMAX SCORE FOR ONTOLOGY "+ontology)
        print(self.keyword_relative_fmax_score)
        # sorted_fmax_rel = sorted(self.keyword_relative_fmax_score[ontology].items(), key=operator.itemgetter(1),
        #  reverse=True)
        sorted_scores = sorted(self.keyword_equal_wts_score.items(), key=operator.itemgetter(1), reverse=True)
        x_ticks = []  # range(0,len(sorted_fmax_rel))
        y1 = []
        y2 = []
        for kwd, val in sorted_scores:
            x_ticks.append(kwd)
            y1.append(val)
            y2.append(self.keyword_relative_fmax_score[kwd])
        # plt.bar(x_index, y1, bar_width, color='#5C89C4', label="Relative FMax")
        # plt.bar(x_index+bar_width, y2, bar_width, color='#5DC687', label="FMax (Equally Weighted)")
        # plt.bar(x_index, y1, bar_width, color='#5C89C4', label="Equal Weights")
        # plt.bar(x_index+bar_width, y2, bar_width, color='#5DC687', label="Weighted by Fmax")
        plt.bar(x_index, y1, bar_width, color='#'+bar1_hex_color, label="Equal Weights")
        plt.bar(x_index+bar_width, y2, bar_width, color='#'+bar2_hex_color, label="Weighted by Fmax")

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
    def csv_ontology_score_results(self):

        # Create CVS Header
        csv_header = list(const.METHODOLOGY_KEYWORDS)
        csv_header.insert(0, "Team")
        csv_header.insert(1, "Taxon ID")
        csv_header.insert(2, "Model Number")
        csv_header.insert(3, "Fmax in "+self.ontology+"_"+"t"+str(self.bnchmrk_type)+"_m"+str(self.bnchmrk_mode))
        csv_header.insert(4, "Equal weights")
        csv_header.insert(5, "Relative Fmax")

        # Create csv input and output file names

        file_name = self.fmax_files_directory + '/' + self.ontology.lower() + "_" + self.taxon_name + '_' + "type" + \
            str(self.bnchmrk_type) + '_' + 'mode' + str(self.bnchmrk_mode) + '_all_fmax_sheet.csv'

        output_file_name = self.csv_output_files_directory + '/' + self.ontology.lower() + "_" + self.taxon_name + "_type"\
            + str(self.bnchmrk_type) + '_' + 'mode' + str(self.bnchmrk_mode) + '_all_fmax_sheet_with_weights.csv'

        # open both input and output csv file
        with open(file_name, newline='') as csvfile_in, open(output_file_name, 'w', newline='') as csvfile_out:

            reader = csv.DictReader(csvfile_in)
            writer = csv.DictWriter(csvfile_out, fieldnames=csv_header)

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
                    if self.taxon_name.lower() == 'all':
                        print("C \'" + self.taxon_name + "\'")

                        if team in self.author_list:
                            for taxon_id in self.author_list[team]:
                                if model in self.author_list[team][str(taxon_id)]:

                                    out_dict.clear()
                                    out_dict['Team'] = team
                                    out_dict['Taxon ID'] = taxon_id
                                    out_dict['Model Number'] = model
                                    out_dict["Fmax in " + self.ontology + "_" + "t" + str(self.bnchmrk_type) + "_m" +
                                             str(self.bnchmrk_mode)] = fmax
                                    out_dict["Equal weights"] = 1.0 / self.total_num_fmax_scores
                                    out_dict['Relative Fmax'] = float(fmax) / self.fmax_sum

                                    for kwd in self.author_list[team][str(taxon_id)][model]:
                                        out_dict[kwd] = 1

                                    writer.writerow(out_dict)
                        else:
                            raise KeyError(team + " not found in author_list for taxonID: " + str(taxon_id))

                    else:
#                        print("D \'" + self.taxon_name + "\'")
#                        print("THIS IS MY AUTHOR LIST")
#                        print(self.author_list)
#                        print("Author TAXON")
#                        print(self.author_list[team][str(self.taxon_id)])

                        if team in self.author_list:
                            if str(self.taxon_id) in self.author_list[team]:
                                if model in self.author_list[team][str(self.taxon_id)]:

                                    out_dict.clear()
                                    out_dict['Team'] = team
                                    out_dict['Taxon ID'] = str(self.taxon_id)
                                    out_dict['Model Number'] = model
                                    out_dict["Fmax in " + self.ontology + "_" + "t" + str(self.bnchmrk_type) + "_m" +
                                             str(self.bnchmrk_mode)] = fmax
                                    out_dict["Equal weights"] = 1.0 / self.total_num_fmax_scores
                                    out_dict['Relative Fmax'] = float(fmax) / self.fmax_sum

                                    for kwd in self.author_list[team][str(self.taxon_id)][model]:
                                        out_dict[kwd] = 1

                                    writer.writerow(out_dict)
                                else:
                                    raise KeyError("Model " + model + " not found in author_list for taxonID: " +
                                                   str(taxon_id) + " Author: " + team)
                            else:
                                raise KeyError("TaxonID: " + str(self.taxon_id) + " not found in author list "
                                                                                  "for author: " + team)
                        else:
                            raise KeyError(team + " not found in author_list for taxonID: " + str(taxon_id))


    def create_taxon_id_dictionary(self):
        # self.taxon_id_counts = {}
        self.taxon_id_counts = Cafa3.get_taxon_id_counts(self.author_list).copy()

    def create_taxon_id_list(self):
        self.taxon_id_list = []

        for taxon_id in self.taxon_id_counts.keys():
            self.taxon_id_list.append(taxon_id)

    def count_kwds(self, model):
        self.model = model
        self.total_keyword_count = self.keyword_template.copy()
        self.model_keyword_count = self.keyword_template.copy()

        print("BEFORE CALC")
        #print(methodology_keyword_count)
        print(self.model_keyword_count)

        # 1. Determine total keyword count by model
        # 2. Determine overall total keyword count
        for author in self.author_list:
            for taxon_id in self.author_list[author]:
                for mdl in self.author_list[author][taxon_id]:
                    for kwrd in self.author_list[author][taxon_id][mdl]:
                        # key_counts['Methodology Keyword Count'] += 1
                        self.total_keyword_count[kwrd] += 1
                        if mdl == model:
                            self.model_keyword_count[kwrd] += 1

        # print("AFTER CALC")
        # print(methodology_keyword_count)
        # print(self.model_keyword_count)

    def plot_raw_keyword_counts(self, width, height, x_font_size, y_font_size, label_font_size, bar_hex_color):
        plt.rc('xtick', labelsize=x_font_size)
        plt.rc('ytick', labelsize=y_font_size)
        plt.rcParams.update({'font.size': label_font_size})

        fig, ax = plt.subplots(figsize=(width, height))
        fig.canvas.set_window_title("Raw Keyword Counts Model " + self.model)

        ax.set_ylabel("Raw Counts")
        ax.set_xlabel("Keyword")

        bar_width = 0.4

        x_index = np.arange(0, len(self.model_keyword_count))

        # Sort in descending order
        sorted_scores = sorted(self.model_keyword_count.items(), key=operator.itemgetter(1), reverse=True)
        x_ticks = []  # range(0,len(sorted_fmax_rel))
        y1 = []
        for kwd, val in sorted_scores:
            x_ticks.append(kwd)
            y1.append(val)
        plt.bar(x_index, y1, bar_width, color='#'+bar_hex_color, label="Equal Weights")

        plt.legend()
        plt.xticks(x_index, x_ticks)
        plt.xticks(rotation=90)
        plt.tight_layout()
        plt.show()

    # Creates two csv files. First csv file outputs all the authors' taxonIDs and models and shows the keywords used
    # The second csv file output the total keywords counts and the total keyword counts by model
    def create_keycount_csv(self):
        csv_header = list(const.METHODOLOGY_KEYWORDS)
        csv_header.insert(0, "Team")
        csv_header.insert(1, "Taxon ID")
        csv_header.insert(2, "Model Number")

        # Create csv file of keywords used for each team's taxonID and model #
        with open(const.CSV_OUTPUT_DIRECTORY+'/team_model_taxonID_keyword.csv', 'w', newline='') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=csv_header)
            writer.writeheader()
            for author in self.author_list:
                for taxonID in self.author_list[author]:
                    for model in self.author_list[author][taxonID]:
                        outline = self.author_list[author][taxonID][model].copy()
                        outline['Team'] = author
                        outline['Taxon ID'] = taxonID
                        outline['Model Number'] = model
                        writer.writerow(outline)  # end of 1st csv file creation

        # Create csv file for total keyword counts and total keywords count by model #
        self.count_kwds('1')
        csv_header = list(const.METHODOLOGY_KEYWORDS)
        csv_header.insert(0, "")
        with open(const.CSV_OUTPUT_DIRECTORY+'/total_keyword_counts.csv', 'w', newline='') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=csv_header)
            writer1 = csv.writer(csvfile, delimiter=' ', quotechar="", quoting=csv.QUOTE_NONE)

            writer1.writerow('Total_Keyword_Count')
            writer.writeheader()
            writer.writerow(self.total_keyword_count)

            csv_header[0] = 'Model Number'
            writer = csv.DictWriter(csvfile, fieldnames=csv_header)
            writer1.writerow('')
            writer1.writerow('Model_Keyword_Count')
            writer.writeheader()
            for i in ('1', '2', '3'):
                self.count_kwds(i)
                row = self.model_keyword_count
                row['Model Number'] = i
                writer.writerow(row)
