import os
import constant
import tarfile


# Directory walk the data folder and place all the CAFA filenames in a Python nested dictionary with structure as
# follows { filename : {'path': filepath} }
def create_file_list(fileList):
    fileCount = 0;
    for dirpath, dirnames, filenames in os.walk(constant.DATA_ROOT_DIRECTORY):
        for file in filenames:
            fileList[file] = {'path': dirpath}
            fileCount += 1
    return fileCount


# Directory walk the a folder and place all the CAFA filenames in a Python List.
def create_file_list_from_tar(fileList):
    fileCount = 0
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
            fileList[filename] = {'path' : path }
            fileCount += 1
            print(tarinfo.name)
            print(filename)
            print(path)
            print(filename,": ",fileList[filename],"\n")

    tarf.close()
    return fileCount
