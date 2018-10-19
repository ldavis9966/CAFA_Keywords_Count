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
        cur_line = fh.readline()
        # if type(cur_line) is bytes:
        #    cur_line = cur_line.decode('utf-8')

        if 'AUTHOR' not in cur_line:
            print("\nError: AUTHOR not found in file.")
            printfile(file_list[file]['path'], file)
            break

        author_in_file = cur_line[len('AUTHOR '):].lower()
        author_in_file = author_in_file.strip()
        author_in_file = author_in_file.replace("\n", "")

        # See if MODEL is in the second line. Extract Model # if so.
        cur_line = fh.readline()
        # if type(cur_line) is bytes:
        #    cur_line = cur_line.decode('utf-8')

        if 'MODEL' not in cur_line:
            print("\nError: MODEL not found in file.")
            printfile(file_list[file]['path'], file)
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
            print("\nError: filename conflicts with file data")
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
            print("\nError: KEYWORDS tag not found in file.")
            printfile(file_list[file]['path'], file)
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

        # Check if the keywords in this file are one of the CAFA accepted keywords
        for kwrd in keywords:
            if kwrd not in const.METHODOLOGY_KEYWORDS:
                print("\nError: keyword", "\""+kwrd+"\"", "not an acceptable methodology keyword in file:")
                printfile(file_list[file]['path'], file)
                break

        # This code block checks if current author is not in the author list. If not present then create a new nested
        # dictionary of both the taxon_id and then another nested dictionary using the current model #.
        if author not in authors_list:
            authors_list[author] = {taxon_id: {model: {}}}
        # Otherwise the author already has an entry. Now check if the current taxon_id is present. If not then create
        # a double nested dictionary for the author using the current taxon_id and model #.
        elif taxon_id not in authors_list[author]:
            authors_list[author][taxon_id] = {model: {}}
        # Otherwise the author already has a taxon_id entry. Now check if the current model has already been
        # added for the current taxon_id. If not, then create a new dictionary for this taxon_id using
        # the current model #.
        elif model not in authors_list[author][taxon_id]:
            authors_list[author][taxon_id][model] = {}
        # If all the prior checks fail, then there must be a double entry of the model # for current taxon_id. So report
        # error. In this case the current Model # will overwrite the prior model # for this taxon_id. Probably should
        # handle this better in the future.
        else:
            print("\nError: duplicate Model # found for Taxon ID")
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


def printfile(path, filename):
    print("File Path + Name:", path+"/"+filename)
    return


# Returns a dictionary where the number of occurences for each taxon_id is of the form
# the {Key: value} = {taxon_id: count}
def get_taxon_id_counts(author_list):
    taxon_id_counts = {}
    for author in author_list:
        for taxon_id in author_list[author]:
            if taxon_id not in taxon_id_counts:
                taxon_id_counts[taxon_id] = 1
            else:
                taxon_id_counts[taxon_id] += 1
    return taxon_id_counts