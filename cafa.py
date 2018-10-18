import file_list as fl
import author_list as al
import taxon
import csv
import operator
import matplotlib.pyplot as plt
import numpy as np
import constant as const
import re


class Cafa:

    def __init__(self):
        self.keyword_template = {
            'sequence alignment': 0, 'sequence-profile alignment': 0, 'profile-profile alignment': 0,
            'phylogeny': 0, 'sequence properties': 0, 'physicochemical properties': 0, 'predicted properties': 0,
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
        self.taxon_id = 0

        self.taxon_id_list = []
        self.taxon_id_counts = {}

        self.set_path_cafa_team_files(const.DATA_ROOT_DIRECTORY)
        self.set_path_zipped_files(const.TAR_DATA_ROOT_DIRECTORY)
        self.set_path_fmax_files(const.FMAX_FILES_DIRECTORY)
        self.set_path_csv_out(const.CSV_OUTPUT_DIRECTORY)
        self.bnchmrk_mode = 1   # Default mode.
        self.bnchmrk_type = 1   # Default type.

        # initialize relative fmax scores dict and scores by equal weighting
        self.keyword_relative_fmax_score = {}
        self.keyword_equal_wts_score = {}
        self.keyword_relative_fmax_score = self.keyword_template.copy()
        self.keyword_equal_wts_score = self.keyword_template.copy()

        self.author_list = {}

        self.fmax_sum = 0
        self.total_num_fmax_scores = 0

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
        self.bnchmrk_type = 1
        self.bnchmrk_mode = 1

        if 'type' in kwargs:
            self.bnchmrk_type = kwargs['type']
        if'mode' in kwargs:
            self.bnchmrk_mode = kwargs['mode']

        # Get either taxonName or taxon_id from argument list
        if'taxonName' in kwargs and 'taxonID' in kwargs:
            print("Both taxonName and taxonID provided as arguments, only one is needed. Ignoring taxonID and using" +
                  " taxonName")
            self.taxon_name = kwargs['taxonName']
            self.taxon_id = taxon.taxon_ID_converter(self.taxon_name)
        elif 'taxonName' in kwargs:
            self.taxon_name = kwargs['taxonName']
            self.taxon_id = taxon.taxon_ID_converter(self.taxon_name)
        elif 'taxonID' in kwargs:
            self.taxon_name = taxon.taxon_name_converter(kwargs['taxonID'])
            self.taxon_id = kwargs['taxonID']
        else:
            raise Exception("'Either taxonName or taxonID needed as an argument.'")

        # taxon_id = taxon.taxon_ID_converter(self.taxon_name)

        # print("Ontology inside tabulator: "+o)
#        if ontology.lower == 'all':
#            file_name = self.fmax_files_directory + '/' + ontology.lower() + "_all_type" + str(bnchmrk_type) + '_' +\
#                    'mode' + str(bnchmrk_mode) + '_all_fmax_sheet.csv'
#        else:
        file_name = self.fmax_files_directory + '/' + ontology.lower() + "_" + self.taxon_name + '_' + "type" +\
            str(self.bnchmrk_type) + '_' + 'mode' + str(self.bnchmrk_mode) + '_all_fmax_sheet.csv'

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
                            if str(taxon_id) in self.author_list[team]:
                                if model in self.author_list[team][str(taxon_id)]:
                                    for kwd in self.author_list[team][str(taxon_id)][model]:
                                        self.keyword_relative_fmax_score[kwd] += float(row['F1-max'])/self.fmax_sum
                                        self.keyword_equal_wts_score[kwd] += 1.0/self.total_num_fmax_scores
                                else:
                                    raise KeyError("Model " + model + " not found in author_list for taxonID: " +
                                                   str(taxon_id) + " Author: " + team)
                            else:
                                raise KeyError("TaxonID: " + str(taxon_id) + "not found in author list for author: " + team)
                        else:
                            raise KeyError(team + " not found in author_list for taxonID: " + str(taxon_id))
                        # print(str(i)+":"+team + " not found in author_list for taxonID: " + str(taxon_id))

    ''' Plots the relative frequency and the relative frequency weighted by fmax   
    '''
    def plot_ontology_fmax_results(self, ontology, **kwargs):

        self.bnchmrk_type = 1
        self.bnchmrk_mode = 1

        ontology = ontology.upper()

        if 'type' in kwargs:
            self.bnchmrk_type = kwargs['type']
        if'mode' in kwargs:
            self.bnchmrk_mode = kwargs['mode']

#        if('taxonName' not in kwargs and 'taxonID' not in kwargs):
#            raise KeyError("Required argument missing. Need either taxonName or taxonID")
        # Get either taxonName or taxonID from argument list
        if 'taxonName' in kwargs and 'taxonID' in kwargs:
            print("Both taxonName and taxonID provided as arguments, only one is needed. Ignoring taxonID and "+\
                  "using taxonName")
            self.taxon_name = kwargs['taxonName']
            self.taxon_id = taxon.taxon_ID_converter(self.taxon_name)
        elif 'taxonName' in kwargs:
            self.taxon_name = kwargs['taxonName']
            self.taxon_id = taxon.taxon_ID_converter(self.taxon_name)
        elif 'taxonID' in kwargs:
            self.taxon_name = taxon.taxon_name_converter(kwargs['taxonID'])
            self.taxon_id = kwargs['taxonID']
        else:
            raise Exception("'Either taxonName or taxonID needed as an argument.'")

        plt.rc('xtick', labelsize=8)
        # fig, ax = plt.subplots(figsize=(10,2))
        # fig = plt.figure(figsize=(10,6))
        fig, ax = plt.subplots(figsize=(8, 6))
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
    def csv_ontology_fmax_results(self, ontology, bnchmrk_type, bnchmrk_mode, taxon_id):

        # Create CVS Header
        csvHeader = list(const.METHODOLOGY_KEYWORDS)
        csvHeader.insert(0, "Team")
        csvHeader.insert(1, "Taxon ID")
        csvHeader.insert(2, "Model Number")
        csvHeader.insert(3, "Fmax in "+ontology+"_"+"t"+str(bnchmrk_type)+"_m"+str(bnchmrk_mode))
        csvHeader.insert(4, "Equal weights")
        csvHeader.insert(5, "Relative Fmax")

        # Create csv input and output file names

        file_name = self.fmax_files_directory + '/' + ontology.lower() + "_" + self.taxon_name + '_' + "type" + \
            str(bnchmrk_type) + '_' + 'mode' + str(bnchmrk_mode) + '_all_fmax_sheet.csv'

        output_file_name = self.csv_output_files_directory + '/' + ontology.lower() + "_" + self.taxon_name + '_' + \
                           "type" + str(bnchmrk_type) + '_' + 'mode' + str(bnchmrk_mode) +\
                           '_all_fmax_sheet_with_weights.csv'

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
                    if self.taxon_name.lower() == 'all':
                        print("C \'" + self.taxon_name + "\'")

                        if team in self.author_list:
                            for(taxon_id) in self.author_list[team]:
                                if model in self.author_list[team][str(taxon_id)]:

                                    out_dict.clear()
                                    out_dict['Team'] = team
                                    out_dict['Taxon ID'] = taxon_id
                                    out_dict['Model Number'] = model
                                    out_dict["Fmax in " + ontology + "_" + "t" + str(bnchmrk_type) + "_m" + str(bnchmrk_mode)] = fmax
                                    out_dict["Equal weights"] = 1.0 / self.total_num_fmax_scores
                                    out_dict['Relative Fmax'] = float(fmax) / self.fmax_sum

                                    for kwd in self.author_list[team][str(taxon_id)][model]:
                                        out_dict[kwd] = 1
                        else:
                            raise KeyError(team + " not found in author_list for taxonID: " + str(taxon_id))

                    else:
                        print("D \'" + self.taxon_name + "\'")
                        if team in self.author_list:
                            if str(taxon_id) in self.author_list[team]:
                                if model in self.author_list[team][str(taxon_id)]:

                                    out_dict.clear()
                                    out_dict['Team'] = team
                                    out_dict['Taxon ID'] = taxon_id
                                    out_dict['Model Number'] = model
                                    out_dict["Fmax in " + ontology + "_" + "t" + str(bnchmrk_type) + "_m" + str(bnchmrk_mode)] = fmax
                                    out_dict["Equal weights"] = 1.0 / self.total_num_fmax_scores
                                    out_dict['Relative Fmax'] = float(fmax) / self.fmax_sum

                                    for kwd in self.author_list[team][str(taxon_id)][model]:
                                        out_dict[kwd] = 1
                                else:
                                    raise KeyError("Model " + model + " not found in author_list for taxonID: " +
                                                   str(taxon_id) + " Author: " + team)
                            else:
                                raise KeyError("TaxonID: " + str(taxon_id) + "not found in author list for author: " + team)
                        else:
                            raise KeyError(team + " not found in author_list for taxonID: " + str(taxon_id))

                    writer.writerow(out_dict)

    def create_taxon_id_dictionary(self):
        # self.taxon_id_counts = {}
        self.taxon_id_counts = al.taxonID_counts(self.author_list).copy()

    def create_taxon_id_list(self):
        self.taxon_id_list = []

        for taxon_id in self.taxon_id_counts.keys():
            self.taxon_id_list.append(taxon_id)

    '''def create_author_list(self, authors_list, file_list):

        files_processed = 0
        for file in file_list:

            # The following code would allow reading the raw submission files directly for a tar file.
            # It was horribly slow taking hours to parse the tar filed and is uncommented for now.
            # May revisit it in future if there is some way to boost performace
            if const.USE_ZIPPED_FILES:
                tarf = tarfile.open(const.TAR_DATA_ROOT_DIRECTORY + "/" + const.TAR_FILE_NAME)
                fh = tarf.extractfile(file_list[file]['path'] + "/" + file)
            else:
                fh = open(file_list[file]['path'] + "/" + file, "r")
            fh = open(file_list[file]['path'] + "/" + file, "r")

            # See if AUTHOR is in the first line. Extract the team name if so.
            curLine = fh.readline()
            if type(curLine) is bytes:
                curLine = curLine.decode('utf-8')

            if 'AUTHOR' not in curLine:
                print("\nError: AUTHOR not found in file.")
                printfile(file_list[file]['path'], file)
                break

            author_in_file = curLine[len('AUTHOR '):].lower()
            author_in_file = author_in_file.strip()
            author_in_file = author_in_file.replace("\n", "")

            # See if MODEL is in the second line. Extract Model # if so.
            curLine = fh.readline()
            if type(curLine) is bytes:
                curLine = curLine.decode('utf-8')

            if 'MODEL' not in curLine:
                print("\nError: MODEL not found in file.")
                printfile(file_list[file]['path'], file)
                break

            model_in_file = curLine[len('MODEL '):]
            model_in_file = model_in_file.replace("\n", "")
            model_in_file = model_in_file.strip(" ")

            fileStrSplit = file.split('_')
            author = fileStrSplit[0].strip().lower()
            model = fileStrSplit[1].strip()
            taxon_id = fileStrSplit[2].split('.')[0].strip(" ")

            # Check if model # and Author (Team Name) are the 1st two items in the file name per CAFA specifications
            if author_in_file.lower() != fileStrSplit[0].lower() or model_in_file != fileStrSplit[1]:
                print("\nError: filename conflicts with file data")
                print("\tFile: " + file_list[file]['path'] + "/" + file)
                # print("\tFilename string: ", file)
                # print("\tFile string parse:", fileStrSplit)
                if author_in_file.lower() != fileStrSplit[0].lower():
                    print("\t\tAuthor in filename does not match author shown in file")
                    print("\t\t\tAuthor (file name):", fileStrSplit[0])
                    print("\t\t\tAuthor (in file):", author_in_file)
                    # Assume the file name author is the correct author name
                    # author = fileStrSplit[0]
                if model_in_file != fileStrSplit[1]:
                    print("\t\tModel # in filename does not match model # shown in file")
                    print("\t\t\tModel (file name):", fileStrSplit[1])
                    print("\t\t\tModel (in file):", model_in_file)
                    # Assume the file name model # is the correct model number
                    # model = fileStrSplit[1]

            curLine = fh.readline()
            if type(curLine) is bytes:
                curLine = curLine.decode('utf-8')

            if 'KEYWORDS' not in curLine:
                print("\nError: KEYWORDS tag not found in file.")
                #printfile(file_list[file]['path'], file)
                break

            # 3 part process: str.split() last part of current line(everything after the KEYWORDS tag) into a list using
            # comma delimiter to extract the keywords. Next remove all characters except upper & lower case
            # a lpha chars, spaces, & hyphens from each keyword in the list. Finally str.strip() any leading 
            # or trailing spaces from each keyword.
            keywords = curLine[len('KEYWORDS '):]
            keywords = keywords.split(",")
            regex = re.compile('[^a-zA-Z -]')
            for i, kwd in enumerate(keywords):
                keywords[i] = regex.sub('', kwd)
                keywords[i] = keywords[i].strip(" ")

            # Check if the keywords in this file are one of the CAFA accepted keywords
            for kwrd in keywords:
                if kwrd not in const.METHODOLOGY_KEYWORDS:
                    print("\nError: keyword", "\"" + kwrd + "\"", "not an acceptable methodology keyword in file:")
                    printfile(file_list[file]['path'], file)
                    break

            # This code block checks if current author is not in the author list. If not present then create a new 
            # nested dictionary of both the taxon_id and then another nested dictionary using the current model #.
            # if author not in authors_list:
            #    authors_list[author] = {taxon_id: {bnchmrk_model: {}}}
            # Otherwise the author already has an entry. Now check if the current taxon_id is present. If not then
            # create a double nested dictionary for the author using the current taxon_id and model #.
            elif taxon_id not in authors_list[author]:
                authors_list[author][taxon_id] = {model: {}}
            # Otherwise the author already has a taxon_id entry. Now check if the current model has already been added
            # for the current taxon_id. If not, then create a new dictionary for this taxon_id using
            # the current model #.
            elif model not in authors_list[author][taxon_id]:
                authors_list[author][taxon_id][model] = {}
            # If all the prior checks fail, then there must be a double entry of the model # for current taxon_id.
            # So report error. In this case the current Model # will overwrite the prior model # for this taxon_id.
            # Probably should handle this better in the future.
            else:
                print("\nError: duplicate Model # found for Taxon ID")
                print("\tTaxonID:", taxon_id)
                print("\tModel #:", model)
                print("\tConflicting file:")
                print("\t\t", file)

            # For every keyword in this current text file, increment the keyword count and store it in the
            # authors_list for this particular author. Also update the total global keyword count
            for kwrd in keywords:
                authors_list[author][taxon_id][model][kwrd] = 1

            # Close current file and increase counter
            fh.close()
            files_processed += 1
        return files_processed


    # Simply helper method that prints a pat and file namme
    def printfile(self, path, filename):
        print("File Path + Name:", path + "/" + filename)


    # Returns a dictionary where the number of each taxon_id is the {Key: value} = {taxon_id: count}
    def taxon_id_counts(self, author_list):
        taxon_id_counts = {}
        for author in author_list:
            for taxon_id in author_list[author]:
                if taxon_id not in taxon_id_counts:
                    taxon_id_counts[taxon_id] = 1
                else:
                    taxon_id_counts[taxon_id] += 1
        return taxon_id_counts'''


'''     # Pre-condition:    
        # 1. Must have called create_author_kwd_taxon_dict() to create the author dictionary
        # 2. Must have called set_path_fmax_files
        # OLD CODE'''
'''        def fmax_kwdscores_by_taxon_id(self, bnchmrk_type, bnchmrk_mode, taxon_id):
            self.keyword_fmax_score_template = {
                'sequence alignment': 0, 'sequence-profile alignment': 0, 'profile-profile alignment': 0,
                'phylogeny': 0,'sequence properties': 0, 'physicochemical properties': 0, 'predicted properties': 0,
                'protein interactions': 0, 'gene expression': 0, 'mass spectrometry': 0, 'genetic interactions': 0,
                'protein structure': 0, 'literature': 0, 'genomic context': 0, 'synteny': 0, 'structure alignment': 0,
                'comparative model': 0, 'predicted protein structure': 0, 'de novo prediction': 0, 
                'machine learning': 0, 'genome environment': 0, 'operon': 0, 'ortholog': 0, 'paralog': 0, 'homolog': 0,
                'hidden Markov model': 0, 'clinical data': 0, 'genetic data': 0, 'natural language processing': 0,
                'other functional information': 0}

            self.keyword_relative_fmax_score = {}
            self.keyword_fmax_score = {}

            for o in const.ONTOLOGY_LIST:
                self.keyword_relative_fmax_score[o] = self.keyword_fmax_score_template.copy()
                self.keyword_fmax_score[o] = self.keyword_fmax_score_template.copy()

            self.taxon_name = taxon.taxon_name_converter(taxon_id)

            self.fmax_sum = {}
            self.total_num_fmax_scores = {}
            for o in const.ONTOLOGY_LIST:
                #print("Ontology inside tabulator: "+o)
                file_name = self.fmax_files_directory + '/' + o.lower() + "_" + self.taxon_name + '_' + "type" + 
                str(bnchmrk_type) + '_' + 'mode' + str(bnchmrk_mode) + '_all_fmax_sheet.csv'

                self.fmax_sum[o] = 0
                self.total_num_fmax_scores[o] = 0
                with open(file_name, newline='') as csvfile:
                    reader = csv.DictReader(csvfile)
                    for i, row in enumerate(reader):
                        if i > 1 and float(row['Coverage']) != 0.0:   # Skip if row is 1 of first 2 rows or coverage = 0
                            self.fmax_sum[o] += float(row['F1-max'])
                            self.total_num_fmax_scores[o] += 1

                with open(file_name, newline='') as csvfile:
                    reader = csv.DictReader(csvfile)
                    for i, row in enumerate(reader):
                        if i > 1 and float(row['Coverage']) != 0.0:   # Skip if row is 1 of first 2 rows or coverage = 0
                            team = row['ID-model'][:-2].lower()
                            model = row['ID-model'][len(row['ID-model'])-1:]

                            if team in self.author_list:
                                if str(taxon_id) in self.author_list[team]:
                                    if model in self.author_list[team][str(taxon_id)]:
                                        for kwd in self.author_list[team][str(taxon_id)][model]:
                                            self.keyword_relative_fmax_score[o][kwd] += 
                                                float(row['F1-max'])/self.fmax_sum[o]
                                            self.keyword_fmax_score[o][kwd] += 1.0/self.total_num_fmax_scores[o]
                                    else:
                                        raise KeyError("Model " + model + " not found in author_list for taxon_id: " +
                                                       str(taxon_id) + " Author: " + team)
                                else:
                                    raise KeyError("taxon_id: " + str(taxon_id) + "not found in author list for author\
                                    : " + team)
                            else:
                                raise KeyError(team + " not found in author_list for taxon_id: " + str(taxon_id))
                                #print(str(i)+":"+team + " not found in author_list for taxon_id: " + str(taxon_id))

            # Debug Code
            #for o in const.ONTOLOGY_LIST:
                #print("Ontology: "+o)
                #print("RELATIVE FMAX")
                #print(self.keyword_relative_fmax_score[o])
                #print("FMAX")
                #print(self.keyword_fmax_score[o])
                sorted_fmax_rel = sorted(self.keyword_relative_fmax_score[o].items(), key=operator.itemgetter(1),
                                    reverse=True)
                #sorted_fmax = sorted(self.keyword_fmax_score[o], key=operator.itemgetter(2), reverse=True)
                #print("SORTED FMAX")
                #print(sorted_fmax)
                #for item in sorted_fmax_rel:
                    #print("%s ---> %.2f" % ( item[0], item[1]))

                #for item in sorted_fmax:
                #    print(item, self.keyword_fmax_score[o][item])
    '''


'''   Plots the relative frequency and the relative frequency weighted by fmax   

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
        print("KEYWORD RELATIVE FMAX SCORE FOR ONTOLOGY "+ontology)
        print(self.keyword_relative_fmax_score[ontology])
        #sorted_fmax_rel = sorted(self.keyword_relative_fmax_score[ontology].items(), key=operator.itemgetter(1),
         reverse=True)
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


'''       Precondition: Must have called self.fmax_kwdscores_by_taxon_id() first.

    def csv_ontology_fmax_results(self, bnchmrk_type, bnchmrk_mode, taxon_id):

        for o in const.ONTOLOGY_LIST:
            # Create CVS Header
            csvHeader = list(const.METHODOLOGY_KEYWORDS)
            csvHeader.insert(0, "Team")
            csvHeader.insert(1, "Taxon ID")
            csvHeader.insert(2, "Model Number")
            csvHeader.insert(3, "Fmax in "+o+"_"+"t"+str(bnchmrk_type)+"_m"+str(bnchmrk_mode))
            csvHeader.insert(4, "Equal weights")
            csvHeader.insert(5, "Relative Fmax")


            # Create csv input and output file names
            file_name = self.fmax_files_directory + '/' + o.lower() + "_" + self.taxon_name + '_' + "type" + 
            str(bnchmrk_type) + '_' + 'mode' + str(bnchmrk_mode) + '_all_fmax_sheet.csv'

            output_file_name = self.csv_output_files_directory + '/' + o.lower() + "_" + self.taxon_name + '_' + "type"
             + str(bnchmrk_type) + '_' + 'mode' + str(bnchmrk_mode) + '_all_fmax_sheet_with_weights.csv'

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
                            if str(taxon_id) in self.author_list[team]:
                                if model in self.author_list[team][str(taxon_id)]:

                                    out_dict.clear()
                                    out_dict['Team'] = team
                                    out_dict['Taxon ID'] = taxon_id
                                    out_dict['Model Number'] = model
                                    out_dict["Fmax in " + o + "_" + "t" + str(bnchmrk_type) + "_m" + 
                                        str(bnchmrk_mode)] = fmax
                                    out_dict["Equal weights"] = 1.0 / self.total_num_fmax_scores[o]
                                    out_dict['Relative Fmax'] = float(fmax) / self.fmax_sum[o]

                                    for kwd in self.author_list[team][str(taxon_id)][model]:
                                        out_dict[kwd] = 1
                                else:
                                    raise KeyError("Model " + model + " not found in author_list for taxon_id: " +
                                                   str(taxon_id) + " Author: " + team)
                            else:
                                raise KeyError("TaxonID: " + str(taxon_id) + "not found in author list for author: " + 
                                team)
                        else:
                            raise KeyError(team + " not found in author_list for taxonID: " + str(taxon_id))

                        writer.writerow(out_dict)'''


'''    def fmax_kwdscores_by_ontology_taxon_list(self, bnchmrk_type, bnchmrk_mode, **kwargs):

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
#            print("taxa: " + k + "\tconverted: " + taxon.taxon_name_converter(k) + "\tunconverted back :" 
#            + taxon.taxon_ID_converter(taxon.taxon_name_converter(k)))

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

#        self.create_taxon_id_dictionary()
#        self.create_taxon_id_list()

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

                #taxon_name = taxon.taxon_name_converter(taxon_id)
                file_name = self.fmax_files_directory + '/' + o.lower() + "_" + taxon_name + '_' + "type" + 
                str(bnchmrk_type) + '_' +\
                            'mode' + str(bnchmrk_mode) + '_all_fmax_sheet.csv'

                with open(file_name, newline='') as csvfile:
                    reader = csv.DictReader(csvfile)
                    for i, row in enumerate(reader):
                        if i > 1 and float(row['Coverage']) != 0.0:   # Skip if row is 1 of first 2 rows or coverage = 0
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
                        if i > 1 and float(row['Coverage']) != 0.0:   # Skip if row is 1 of first 2 rows or coverage = 0
                            team = row['ID-model'][:-2].lower()
                            model = row['ID-model'][len(row['ID-model'])-1:]

                            if team in self.author_list:
                                if str(taxon_id) in self.author_list[team]:
                                    if model in self.author_list[team][str(taxon_id)]:
                                        for kwd in self.author_list[team][str(taxon_id)][model]:
                                            keyword_relative_fmax_score[o][kwd] += float(row['F1-max']) / fmax_sum[o]
                                            keyword_equal_wts_score[o][kwd] += 1.0 / total_num_fmax_scores[o]
                                            keyword_relative_fmax_score['ALL'][kwd] += 
                                                    float(row['F1-max']) / fmax_sum_all_ontologies
                                            keyword_equal_wts_score['ALL'][kwd] += 
                                                    1.0 / total_num_fmax_scores_all_ontologies

                                    else:
                                        raise KeyError("Model " + model + " not found in author_list for taxonID: " +
                                                       str(taxon_id) + " Author: " + team)
                                else:
                                    raise KeyError("TaxonID: " + str(taxon_id) + "not found in author list for author: "
                                        + team)
                            else:
                                raise KeyError(team + " not found in author_list for taxonID: " + str(taxon_id))
                                #print(str(i)+":"+team + " not found in author_list for taxonID: " + str(taxon_id))


        # Debug Code
        #for o in const.ONTOLOGY_LIST:
            #print("Ontology: "+o)
            #print("RELATIVE FMAX")
            #print(self.keyword_relative_fmax_score[o])
            #print("FMAX")
            #print(self.keyword_fmax_score[o])
            sorted_fmax_rel = sorted(self.keyword_relative_fmax_score[o].items(), key=operator.itemgetter(1), 
                reverse=True)
            #sorted_fmax = sorted(self.keyword_fmax_score[o], key=operator.itemgetter(2), reverse=True)
            #print("SORTED FMAX")
            #print(sorted_fmax)
            #for item in sorted_fmax_rel:
                #print("%s ---> %.2f" % ( item[0], item[1]))

            #for item in sorted_fmax:
            #    print(item, self.keyword_fmax_score[o][item]) 
'''