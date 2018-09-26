import re
import tarfile
import constant as const



def create_author_list(authors_list, file_list):

    files_processed = 0

    for file in file_list:
        if const.USE_ZIPPED_FILES:
            tarf = tarfile.open(const.TAR_DATA_ROOT_DIRECTORY + "/" + const.TAR_FILE_NAME)
            fh = tarf.extractfile(file_list[file]['path']+"/"+file)
        else:
            fh = open(file_list[file]['path']+"/"+file, "r")

        # See if AUTHOR is in the first line. Extract the team name if so.
        curLine = fh.readline()
        if type(curLine) is bytes:
            curLine = curLine.decode('utf-8')

        if 'AUTHOR' not in curLine:
            print("\nError: AUTHOR not found in file.")
            printfile(file_list[file]['path'], file)
            break

        author = curLine[len('AUTHOR '):]
        author = author.replace("\n", "")

        # See if MODEL is in the second line. Extract Model # if so.
        curLine = fh.readline()
        if type(curLine) is bytes:
            curLine = curLine.decode('utf-8')

        if 'MODEL' not in curLine:
            print("\nError: MODEL not found in file.")
            printfile(file_list[file]['path'], file)
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
            print("\t\t\t"+file_list[file]['path'], file)

        curLine = fh.readline()
        if type(curLine) is bytes:
            curLine = curLine.decode('utf-8')

        if 'KEYWORDS' not in curLine:
            print("\nError: KEYWORDS tag not found in file.")
            printfile(file_list[file]['path'], file)
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
            if kwrd not in const.METHODOLOGY_KEYWORDS:
                print("\nError: keyword", "\""+kwrd+"\"", "not an acceptable methodology keyword in file:")
                printfile(file_list[file]['path'], file)
                break

        # This code block checks if current author is not in the author list. If not present then create a new nested
        # dictionary of both the taxonID and then another nested dictionary using the current model #.
        if author not in authors_list:
            authors_list[author] = {taxonID: {model: {}}}
        # Otherwise the author already has an entry. Now check if the current taxonID is present. If not then create
        # a double nested dictionary for the author using the current taxonID and model #.
        elif taxonID not in authors_list[author]:
            authors_list[author][taxonID] = {model: {}}
        # Otherwise the author already has a taxonID entry. Now check if the current model has already been added for the
        # current taxonID. If not, then create a new dictionary for this taxonID using the current model #.
        elif model not in authors_list[author][taxonID]:
            authors_list[author][taxonID][model] = {}
        # If all the prior checks fail, then there must be a double entry of the model # for current taxonID. So report
        # error. In this case the current Model # will overwrite the prior model # for this taxonID. Probably should
        # handle this better in the future.
        else:
            print("\nError: duplicate Model # found for Taxon ID")
            print("\tTaxonID:", taxonID)
            print("\tModel #:", model)
            print("\tConflicting file:")
            print("\t\t", file)

        # For every keyword in this current text file, increment the keyword count and store it in the authors_list for this
        # particular author.
        # Also update the total global keyword count
        for kwrd in keywords:
            authors_list[author][taxonID][model][kwrd] = 1

        # Close current file and increase counter
        fh.close()
        files_processed += 1
    return files_processed


def printfile(path, filename):
    print("File Path + Name:", path+"/"+filename)
    return


