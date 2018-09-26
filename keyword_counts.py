import constant

# Note: this dictionary must contain the same exact keywords as in constant.METHODOLOGY_KEYWORDS
methodology_keyword_count = {'sequence alignment': 0,
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
                             'other functional information': 0}

model_methodology_keyword_count = {"1": {}, "2": {}, "3": {}}

def count_kwds(authors_list):
    # Build empty keyword count dictionary by model number
    # When populated has the format: {'model#': {'keyword1': count1, 'keyword2': count2, etc }, etc. }
    # Example below:
    # { '1':{'ortholog': 234, 'sequence-profile alignment':34, etc.}, '2': {'machine learning': 42, 'operon': 12} }
    for kwd in constant.METHODOLOGY_KEYWORDS:
        for i in ['1', '2', '3']:
            model_methodology_keyword_count[i][kwd] = 0

    # 1. Determine total keyword count by model #
    # 2. Determine overall total keyword count
    for author in authors_list:
        for taxonID in authors_list[author]:
            for model in authors_list[author][taxonID]:
                for kwrd in authors_list[author][taxonID][model]:
                    # key_counts['Methodology Keyword Count'] += 1
                    methodology_keyword_count[kwrd] += 1
                    model_methodology_keyword_count[model][kwrd] += 1

