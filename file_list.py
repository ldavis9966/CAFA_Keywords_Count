import os
import constant
import tarfile
import status_bar as sb

# Directory walk the data folders (from either a zipped file or unzipped file) and place all the CAFA filenames in a
# Python nested dictionary with structure as follows { filename : {'path': filepath} }.
def create_file_list(fileList):
    file_count = 0

    if constant.USE_ZIPPED_FILES: # use zipped files
        tarf = tarfile.open(constant.TAR_DATA_ROOT_DIRECTORY + "/" + constant.TAR_FILE_NAME)
        #    tarf = tarfile.open("data.tar.gz")

        for tarinfo in tarf:
            if tarinfo.isreg():
                tmp = tarinfo.name.split("/")
                filename = tmp[len(tmp) - 1]
                del tmp[-1]
                path = ""
                for dir in tmp:
                    path = path + dir + "/"
                path = path.strip("/")
                fileList[filename] = {'path': path}
                file_count += 1

                sb.print(2740, 20, file_count)

#                print(tarinfo.name)
#                print(filename)
#                print(path)
#                print(filename, ": ", fileList[filename], "\n")
        tarf.close()
    else: # Use unzipped files
        for dirpath, dirnames, filenames in os.walk(constant.DATA_ROOT_DIRECTORY):
            for file in filenames:
                fileList[file] = {'path': dirpath}
                file_count += 1

    return file_count