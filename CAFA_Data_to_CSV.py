
import os
import re
import csv
import constant
import file_list as fl
import tarfile

'''
To do list:
1. Refine the error handling some more
2. add functions
3. read from zip file instead of unpacking
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
This code relies on the follow directory structure to function correctly

data
   |_raw_submission
                  |_TEAM NAME(S)
                              |_ Model_#
                                       |_Prediction Files. Filename format should be as follows:
                                                                (TeamID)_(Model #)_(Taxon ID).txt


NOTE: For this program to execute, the DATA_ROOT_DIRECTORY path must be set in the file constant.py. This path must 
        lead to the folder with the TEAM NAMES directory as shown above.  
'''

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


def printfile(path, filename):
    print("File Path + Name:", path+"/"+filename)
    return


fileCount = 0
fileList = {}
authorsList = {}

# Directory walk the data folder and place all the CAFA filenames in a Python List.
#fileCount = fl.create_file_list(fileList)
fileCount = fl.create_file_list_from_tar(fileList)



filesProcessed = 0
for file in fileList:
#    fh = open(fileList[file]['path']+"/"+file, "r")
    tarf = tarfile.open(constant.TAR_DATA_ROOT_DIRECTORY + "/" + constant.TAR_FILE_NAME)
    fh = tarf.extractfile(fileList[file]['path']+"/"+file)




    # See if AUTHOR is in the first line. Extract the team name if so.
    curLine = fh.readline()
    if type(curLine) is bytes:
        curLine = curLine.decode('utf-8')

    if 'AUTHOR' not in curLine:
        print("\nError: AUTHOR not found in file.")
        printfile(fileList[file]['path'], file)
        break

    author = curLine[len('AUTHOR '):]
    author = author.replace("\n", "")

    # See if MODEL is in the second line. Extract Model # if so.
    curLine = fh.readline()
    if type(curLine) is bytes:
        curLine = curLine.decode('utf-8')

    if 'MODEL' not in curLine:
        print("\nError: MODEL not found in file.")
        printfile(fileList[file]['path'], file)
        break

    model = curLine[len('MODEL '):]
    model = model.replace("\n", "")
    model = model.strip(" ")

    fileStrSplit = file.split('_')
    taxonID = fileStrSplit[2].split('.')[0].strip(" ")

    # Check if model # and Author (Team Name) are the 1st two items in the file name per CAFA specifications
    if author not in fileStrSplit and model not in fileStrSplit:
        print("\nError: filename conflicts with file data")
        print("\tFilename string: ", file)
        print("\tFile string parse:", fileStrSplit)
        if author != fileStrSplit[0]:
            print("\t\tAuthor in filename does not match author shown in file")
            print("\t\t\tAuthor (file name):", fileStrSplit[0])
            print("\t\t\tAuthor (in file):", author)
            # Assume the file name author is the correct author name
            author = fileStrSplit[0]
        if model != fileStrSplit[1]:
            print("\t\tModel # in filename does not match model # shown in file")
            print("\t\t\tModel (file name):", fileStrSplit[1])
            print("\t\t\tModel (in file):", model)
            # Assume the file name model # is the correct model number
            model = fileStrSplit[1]
        print("\t\t\t"+fileList[file]['path'], file)

    curLine = fh.readline()
    if type(curLine) is bytes:
        curLine = curLine.decode('utf-8')

    if 'KEYWORDS' not in curLine:
        print("\nError: KEYWORDS tag not found in file.")
        printfile(fileList[file]['path'], file)
        break

    # 3 part process: str.split() last part of current line(everything after the KEYWORDS tag) into a list using
    # comma delimiter to extract the keywords. Next remove all characters except upper & lower case alpha chars, spaces,
    # & hyphens from each keyword in the list. Finally str.strip() any leading or trailing spaces from each keyword.
    keywords = curLine[len('KEYWORDS '):]
    keywords = keywords.split(",")
    regex = re.compile('[^a-zA-Z -]')
    for i, kwd in enumerate(keywords):
        keywords[i] = regex.sub('', kwd)
        keywords[i] = keywords[i].strip(" ")

    # Check if the keywords in this file are one of the CAFA accepted keywords
    for kwrd in keywords:
        if kwrd not in constant.METHODOLOGY_KEYWORDS:
            print("\nError: keyword", "\""+kwrd+"\"", "not an acceptable methodology keyword in file:")
            printfile(fileList[file]['path'], file)
            break

    # This code block checks if current author is not in the author list. If not present then create a new nested
    # dictionary of both the taxonID and then another nested dictionary using the current model #.
    if author not in authorsList:
        authorsList[author] = {taxonID: {model: {}}}
    # Otherwise the author already has an entry. Now check if the current taxonID is present. If not then create
    # a double nested dictionary for the author using the current taxonID and model #.
    elif taxonID not in authorsList[author]:
        authorsList[author][taxonID] = {model: {}}
    # Otherwise the author already has a taxonID entry. Now check if the current model has already been added for the
    # current taxonID. If not, then create a new dictionary for this taxonID using the current model #.
    elif model not in authorsList[author][taxonID]:
        authorsList[author][taxonID][model] = {}
    # If all the prior checks fail, then there must be a double entry of the model # for current taxonID. So report
    # error. In this case the current Model # will overwrite the prior model # for this taxonID. Probably should
    # handle this better in the future.
    else:
        print("\nError: duplicate Model # found for Taxon ID")
        print("\tTaxonID:", taxonID)
        print("\tModel #:", model)
        print("\tConflicting file:")
        print("\t\t", file)

    # For every keyword in this current text file, increment the keyword count and store it in the authorsList for this
    # particular author.
    # Also update the total global keyword count
    for kwrd in keywords:
        authorsList[author][taxonID][model][kwrd] = 1

    # Close current file and increase counter
    fh.close()
    filesProcessed += 1


# Create csv file header row
csvHeader = list(constant.METHODOLOGY_KEYWORDS)
csvHeader.insert(0, "Team")
csvHeader.insert(1, "Taxon ID")
csvHeader.insert(2, "Model Number")

# Create csv file of keywords used for each team's taxonID and model #
with open('team_model_taxonID_keyword.csv', 'w', newline='') as csvfile:
    writer = csv.DictWriter(csvfile, fieldnames=csvHeader)
    writer.writeheader()
    for author in authorsList:
        for taxonID in authorsList[author]:
            for model in authorsList[author][taxonID]:
                outline = authorsList[author][taxonID][model].copy()
                outline['Team'] = author
                outline['Taxon ID'] = taxonID
                outline['Model Number'] = model
                writer.writerow(outline)    # end of 1st csv file creation


# Build empty keyword count dictionary by model number
model_methodology_keyword_count = {"1": {}, "2": {}, "3": {}}
for kwd in constant.METHODOLOGY_KEYWORDS:
    for i in ['1', '2', '3']:
        model_methodology_keyword_count[i][kwd] = 0

# 1. Determine total keyword count by model #
# 2. Determine overall total keyword count
for author in authorsList:
    for taxonID in authorsList[author]:
        for model in authorsList[author][taxonID]:
            for kwrd in authorsList[author][taxonID][model]:
                methodology_keyword_count[kwrd] += 1
                model_methodology_keyword_count[model][kwrd] += 1

# Create csv file for total keyword counts and total keywords count by model #
csvHeader = list(constant.METHODOLOGY_KEYWORDS)
csvHeader.insert(0, "")
with open('total_keyword_counts.csv', 'w', newline='') as csvfile:
    writer = csv.DictWriter(csvfile, fieldnames=csvHeader)
    writer1 = csv.writer(csvfile, delimiter=' ', quotechar="", quoting=csv.QUOTE_NONE)

    writer1.writerow('Total_Keyword_Count')
    writer.writeheader()
    writer.writerow(methodology_keyword_count)

    csvHeader[0] = 'Model Number'
    writer = csv.DictWriter(csvfile, fieldnames=csvHeader)
    writer1.writerow('')
    writer1.writerow('Model_Keyword_Count')
    writer.writeheader()
    for i in ('1', '2', '3'):
        row = model_methodology_keyword_count[i]
        row['Model Number'] = i
        writer.writerow(row)


# Output data
print("\n\n\nAUTHORS LIST DICTIONARY")
print(authorsList)
print("\nMETHODOLOGY KEYWORD COUNT DICTIONARY")
print(methodology_keyword_count)
print("\nMODEL METHODOLOGY KEYWORD COUNT DICTIONARY")
print(model_methodology_keyword_count)
print("\nNumber of files discovered:", fileCount)
print("\nNumber of files processed", filesProcessed)
