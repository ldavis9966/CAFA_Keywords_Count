import os
import constant
import tarfile
import status_bar as sb


# Directory walk the data folders (from either a zipped file or unzipped file) and place all the CAFA filenames in a
# Python nested dictionary with structure as follows { filename : {'path': filepath} }.
def create_file_list(file_list):
    file_count = 0

    if constant.USE_ZIPPED_FILES:  # use zipped files
        tarf = tarfile.open(constant.TAR_DATA_ROOT_DIRECTORY + "/" + constant.TAR_FILE_NAME)
        #    tarf = tarfile.open("data.tar.gz")

        for tarinfo in tarf:
            if tarinfo.isreg():                 # Check tarinfo is a file, and not a directory
                tmp = tarinfo.name.split("/")   # Split out the path+filename strings into a tuple using '/'
                filename = tmp[len(tmp) - 1]    # Get filename from last tuple
                del tmp[-1]                     # last entry from tuple which is the filename.
                path = ""                       # now rebuild the path from the tuple.
                for directory in tmp:
                    path = path + directory + "/"
                path = path.strip("/")
                file_list[filename] = {'path': path}
                file_count += 1

                sb.print(2740, 20, file_count)

#                print(tarinfo.name)
#                print(filename)
#                print(path)
#                print(filename, ": ", file_list[filename], "\n")
        tarf.close()
    else:  # Use unzipped files
        for dirpath, dirnames, filenames in os.walk(constant.DATA_ROOT_DIRECTORY):
            for file in filenames:
                file_list[file] = {'path': dirpath}
                file_count += 1

    return file_count
