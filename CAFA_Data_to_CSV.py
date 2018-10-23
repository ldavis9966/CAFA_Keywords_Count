import file_list as fl
import author_list as al
import CSV_files_generator as cfg
import keyword_counts as kc
'''
To do list:
1. Refine the error handling some more
'''

'''
This program traverses all the CAFA files within a data folder to create several csv files. It also performs several
checks including the following:

1.  Ensure that a filename matches its internal data. Per CAFA specs, a text filename should be formatted to include the
    TeamID, the model number (1,2 or 3), and the taxonID (in that order), with underscores separating each. For 
    example: GOBBI_2_99287.txt Within each file, the first 2 lines of the file contain the Author and the Model number.
    The program parses the file name and compares it to the data inside the file to make sure the model and
    author(TeamID) are the same. 
2.  The third line of a file includes a list of acceptable keywords for the method of protein prediction. The program
    ensures that the file's listed keywords are one of the acceptable keywords.    

The program walks the directory of data and stores each filename in a dictionary. Then the program traverses that
dictionary to create a nested dictionary for each team's submission. The directory structure is constructed as

    {TeamID:{TaxonID:{Model#:{Keyword:1}}}} where only keywords present are given a key/value pair of 1  

Once the nested dictionary is built, it then traverses the dictionary to output to a csv file the keywords used for each
team's submission where each row of the csv is of the format:
    Team, Taxon ID, Model Number, sequence alignment, sequence-profile alignment, profile-profile alignment, (etc.) 

Using the nested dictionary, other outputted data formats can easily be created for different csv files. 

Preconditions;
1. The program can either read in unzipped files or zipped files. The root directory for unzipped files must be
    specified in the constant.py file using the DATA_ROOT_DIRECTORY var.
2. For the zipped files (tar.gz), the zip file name and root directory must be specified in contant.py using the
    TAR_FILE_NAME and TAR_DATA_ROOT_DIRECTORY vars respectively.
3. To use zipped files instead of unzipped, set the flag USE_ZIPPED_FILES to true in constant.py
4. This code relies on the follow directory structure to function correctly:

Directory Structure (for unzippled files)
data
   |_raw_submission
                  |_TEAM NAME(S)
                              |_ Model_#
                                       |_Prediction Files. Filename format should be as follows:
                                                                (TeamID)_(Model #)_(Taxon ID).txt


'''


fileList = {}
authorsList = {}

fileCount = fl.create_file_list(fileList)
filesProcessed = al.create_author_list(authorsList, fileList)

#cfg.create_keycount_csv(authorsList)

#kc.count_kwds(authorsList)
#kc.plot_raw_keyword_counts("none", kc.methodology_keyword_count, taxonID=9606)
print(kc.model_methodology_keyword_count)
#kc.plot_raw_keyword_counts("none", kc.model_methodology_keyword_count['1'], taxonID=9606)
#kc.plot_raw_keyword_counts("none", kc.model_methodology_keyword_count['2'], taxonID=9606)
#kc.plot_raw_keyword_counts("none", kc.model_methodology_keyword_count['3'], taxonID=9606)
key_cnts_for_model = kc.count_kwds(authorsList, '1')
kc.plot_raw_keyword_counts(key_cnts_for_model, '1', 8, 6, 8, 8, 12, 5C89C4)
key_cnts_for_model = kc.count_kwds(authorsList, '2')
kc.plot_raw_keyword_counts(key_cnts_for_model, '2', 8, 6, 8, 8, 12, 5C89C4)
key_cnts_for_model = kc.count_kwds(authorsList, '3')
kc.plot_raw_keyword_counts(key_cnts_for_model, '3', 8, 6, 8, 8, 12, ff0000)
print(kc.model_methodology_keyword_count)

# Output data
#print("\n\n\nAUTHORS LIST DICTIONARY")
#print(authorsList)
#print("\nMETHODOLOGY KEYWORD COUNT DICTIONARY")
#print(kc.methodology_keyword_count)
#print("\nMODEL METHODOLOGY KEYWORD COUNT DICTIONARY")
#print(kc.model_methodology_keyword_count)
#print("\nNumber of files discovered:", fileCount)
#print("\nNumber of files processed", filesProcessed)


#results{}

print("Before")
#for author in authorsList:
#    for taxonID in authorsList[author]:
#        if author == 'casadiobolognabiocomputinggroup' and taxonID == '9606':
#            print(author,taxonID, authorsList[author][taxonID].keys())
#            print(author,taxonID, authorsList[author][taxonID])
print("After")


