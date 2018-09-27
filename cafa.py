import file_list as fl
import author_list as al
import taxon
import csv

class cafa:

    #cafa
    #DATA_ROOT_DIRECTORY = 'data/raw_submission'
    #CSV_OUTPUT_DIRECTORY = 'csv_files'
    #TAR_DATA_ROOT_DIRECTORY = 'data'
    #TAR_FILE_NAME = 'raw_submission.tar.gz'

    def set_path_cafa_team_files(self, path):
        self.data_root_directory = path

    def set_path_zipped_files(self, path):
        self.zip_data_directory = path

    def set_path_fmax_files(self, path):
        self.fmax_files_directory = path

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
    def fmax_kwdscores_by_taxonID(self, ontology, type, mode, taxonID):
        self.keyword_fmax_score = {
            'sequence alignment': 0, 'sequence-profile alignment': 0, 'profile-profile alignment': 0,
            'phylogeny': 0,'sequence properties': 0, 'physicochemical properties': 0, 'predicted properties': 0,
            'protein interactions': 0, 'gene expression': 0, 'mass spectrometry': 0, 'genetic interactions': 0,
            'protein structure': 0, 'literature': 0, 'genomic context': 0, 'synteny': 0, 'structure alignment': 0,
            'comparative model': 0, 'predicted protein structure': 0, 'de novo prediction': 0, 'machine learning': 0,
            'genome environment': 0, 'operon': 0, 'ortholog': 0, 'paralog': 0, 'homolog': 0, 'hidden Markov model': 0,
            'clinical data': 0, 'genetic data': 0, 'natural language processing': 0, 'other functional information': 0}

        ontology = ontology.upper().lower()
        if ( (ontology != 'mfo') and (ontology != 'bpo') and (ontology != 'cco') ):
            raise ValueError("Ontology not recognized. Must be BPO, CCO, or MFO.")

        taxon_name = taxon.taxon_name_converter(taxonID)
        file_name = self.fmax_files_directory + '/' + ontology + "_" + taxon_name + '_' + "type" + str(type) + '_' +\
                    'mode' + str(mode) + '_all_fmax_sheet.csv'


        with open(file_name, newline='') as csvfile:
            reader = csv.DictReader(csvfile)
            for i, row in enumerate(reader):
                if i > 1:   # first 2 rows are non team data
                    team = row['ID-model'][:-2]
                    model = row['ID-model'][len(row['ID-model'])-1:]
                    #print(row['ID-model'], team, model)

        located_dict = self.author_list[team][]
        for author in self.author_list:
            if author == team:
                for dict_taxonID in self.author_list[author]
                    if dict_taxonID == taxonID:
                        for dict_model in self.author_list[author][]
            for taxonID in self.author_list[author]
                for model in self.author_list[author]



        print('ontology:', ontology)
        print('taxonID', taxonID)
        print(file_name)

