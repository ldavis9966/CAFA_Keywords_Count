import constant
import taxon
import matplotlib.pyplot as plt
import numpy as np
import operator

# Note: this dictionary must contain the same exact keywords as in constant.METHODOLOGY_KEYWORDS
'''methodology_keyword_count = {'sequence alignment': 0,
                             'sequence-profile alignment': 0,
                             'profile-profile alignment': 0,
                             'phylogeny': 0,
                             'sequence properties': 0,
                             'physicochemical properties': 0,
                             'predicted properties': 0,
                             'protein interactions': 0,
                             'gene expression': 0,
                             'mass spectrometry': 0,
                             'genetic interactions': 0,
                             'protein structure': 0,
                             'literature': 0,
                             'genomic context': 0,
                             'synteny': 0,
                             'structure alignment': 0,
                             'comparative model': 0,
                             'predicted protein structure': 0,
                             'de novo prediction': 0,
                             'machine learning': 0,
                             'genome environment': 0,
                             'operon': 0,
                             'ortholog': 0,
                             'paralog': 0,
                             'homolog': 0,
                             'hidden Markov model': 0,
                             'clinical data': 0,
                             'genetic data': 0,
                             'natural language processing': 0,
                             'other functional information': 0}'''

keyword_template = {
            'sequence alignment': 0, 'sequence-profile alignment': 0, 'profile-profile alignment': 0,
            'phylogeny': 0, 'sequence properties': 0, 'physicochemical properties': 0, 'predicted properties': 0,
            'protein interactions': 0, 'gene expression': 0, 'mass spectrometry': 0, 'genetic interactions': 0,
            'protein structure': 0, 'literature': 0, 'genomic context': 0, 'synteny': 0, 'structure alignment': 0,
            'comparative model': 0, 'predicted protein structure': 0, 'de novo prediction': 0, 'machine learning': 0,
            'genome environment': 0, 'operon': 0, 'ortholog': 0, 'paralog': 0, 'homolog': 0, 'hidden Markov model': 0,
            'clinical data': 0, 'genetic data': 0, 'natural language processing': 0, 'other functional information': 0}

model_methodology_keyword_count = {"1": {}, "2": {}, "3": {}}

methodology_keyword_count = keyword_template.copy()


def count_kwds(authors_list):

    for model in ['1', '2', '3']:
        model_methodology_keyword_count[model] = keyword_template.copy()

    #print(model_methodology_keyword_count['1'])

    print("BEFORE CALC")
    #print(methodology_keyword_count)
    print(model_methodology_keyword_count)


    # Build empty keyword count dictionary by model number
    # When populated has the format: {'model#': {'keyword1': count1, 'keyword2': count2, etc }, etc. }
    # Example below:
    # { '1':{'ortholog': 234, 'sequence-profile alignment':34, etc.}, '2': {'machine learning': 42, 'operon': 12} }
#    for kwd in constant.METHODOLOGY_KEYWORDS:
#        for model in ['1', '2', '3']:
#            model_methodology_keyword_count[model][kwd] = 0

    # 1. Determine total keyword count by model #
    # 2. Determine overall total keyword count
    for author in authors_list:
        for taxonID in authors_list[author]:
            for model in authors_list[author][taxonID]:
                for kwrd in authors_list[author][taxonID][model]:
                    # key_counts['Methodology Keyword Count'] += 1
                    methodology_keyword_count[kwrd] += 1
                    model_methodology_keyword_count[model][kwrd] += 1

    print("AFTER CALC")
    #print(methodology_keyword_count)
    print(model_methodology_keyword_count)

def plot_raw_keyword_counts(ontology, keyword_dict, **kwargs):

#        self.bnchmrk_type = 1
#        self.bnchmrk_mode = 1

    ontology = ontology.upper()

#        if 'type' in kwargs:
#            self.bnchmrk_type = kwargs['type']
#        if'mode' in kwargs:
#            self.bnchmrk_mode = kwargs['mode']

#        if('taxonName' not in kwargs and 'taxonID' not in kwargs):
#            raise KeyError("Required argument missing. Need either taxonName or taxonID")
    # Get either taxonName or taxonID from argument list
    '''if 'taxonName' in kwargs and 'taxonID' in kwargs:
        print("Both taxonName and taxonID provided as arguments, only one is needed. Ignoring taxonID and " +
              "using taxonName")
        self.taxon_name = kwargs['taxonName']
        self.taxon_id = taxon.taxon_name_to_id(self.taxon_name)
    elif 'taxonName' in kwargs:
        self.taxon_name = kwargs['taxonName']
        self.taxon_id = taxon.taxon_name_to_id(self.taxon_name)
    elif 'taxonID' in kwargs:
        self.taxon_name = taxon.taxon_id_to_name(kwargs['taxonID'])
        self.taxon_id = kwargs['taxonID']
    else:
        raise Exception("'Either taxonName or taxonID needed as an argument.'")'''

    plt.rc('xtick', labelsize=8)
    # fig, ax = plt.subplots(figsize=(10,2))
    # fig = plt.figure(figsize=(10,6))
    fig, ax = plt.subplots(figsize=(8, 6))
#        fig.suptitle("Ontology: "+ontology+"   Taxon: "+ self.taxon_name, fontsize=14)
#        fig.canvas.set_window_title("Ontology: " + ontology + "\tTaxon: " + self.taxon_name+"\tType: " +
#                                    str(self.bnchmrk_type) + "\tMode: " + str(self.bnchmrk_mode))
    fig.canvas.set_window_title("Raw Keyword Counts for All Models")

# ax = fig.add_axes([.1, .25, .8, .7])
    # ax.set_title(ontology + " Cumulative Relative Fmax Scores by Keyword for " + self.taxon_name + " Taxon")
    ax.set_ylabel("Raw Counts")
    ax.set_xlabel("Keyword")

    bar_width = 0.4

    #x_index = np.arange(0, len(self.keyword_relative_fmax_score))
    x_index = np.arange(0, len(keyword_dict))

# Sorts the keywords of the Dictionary, self.keyword_relative_fmax_score, by the values of its keywords
    ####print("KEYWORD RELATIVE FMAX SCORE FOR ONTOLOGY "+ontology)
    ####print(self.keyword_relative_fmax_score)
    # sorted_fmax_rel = sorted(self.keyword_relative_fmax_score[ontology].items(), key=operator.itemgetter(1),
    #  reverse=True)
    sorted_scores = sorted(keyword_dict.items(), key=operator.itemgetter(1), reverse=True)
    x_ticks = []  # range(0,len(sorted_fmax_rel))
    y1 = []
    #y2 = []
    for kwd, val in sorted_scores:
        x_ticks.append(kwd)
        y1.append(val)
        #y2.append(self.keyword_relative_fmax_score[kwd])
    # plt.bar(x_index, y1, bar_width, color='#5C89C4', label="Relative FMax")
    # plt.bar(x_index+bar_width, y2, bar_width, color='#5DC687', label="FMax (Equally Weighted)")
    plt.bar(x_index, y1, bar_width, color='#5C89C4', label="Equal Weights")
    #plt.bar(x_index+bar_width, y2, bar_width, color='#5DC687', label="Weighted by Fmax")

    index = np.arange(4)
    print(index)

    plt.legend()
    plt.xticks(x_index+bar_width/2, x_ticks)
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.show()
    print(x_index)
    print(x_index+bar_width)
