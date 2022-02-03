"""
This script will populate a text file with all ROOT files available in a given CRAB output directory.

# python write_input_file.py --dir /store/user/srosenzw/path/to/files/######_######/0000/
# python write_input_file.py --dir /store/group/lpchbb/srosenzw/path/to/files/######_######/0000/
"""

from argparse import ArgumentParser
import os
import re
import subprocess
# import sys

print(".. parsing argument line")
parser = ArgumentParser(description='Command line parser of model options and tags')
parser.add_argument('--year', dest='year', help='conditions of which year', required=True)
parser.add_argument('--dir', dest='dir', help='CRAB output directory', required=True)

args = parser.parse_args()

dirName = args.dir

# find the desired file name, NMSSM_XYH_YToHH_6b_MX_###_MY_###
# Assumption: crab job was saved in the format
# srosenzw_NMSSM_XYH_YToHH_6b_MX_###_MY_###_sl7_nano_100k
start = re.search('NMSSM*',dirName).start()
end = re.search('_sl7',dirName).start()
textfile = f"{dirName[start:end]}.txt"

def exists_on_eos(lfn):
    """ check if lfn (starting with /store/group) exists """
    retcode = os.system('eos root://cmseos.fnal.gov ls -s %s > /dev/null 2>&1' % lfn)
    # print "THE FOLDER", lfn, "RETURNED CODE", retcode
    return True if retcode == 0 else False

outputName = f"input/PrivateMC_{args.year}/{textfile}"
with open(outputName, "w") as f:
    print(f".. writing to file: {outputName}")
    if (exists_on_eos(dirName)):
        output = subprocess.run(["eos", "root://cmseos.fnal.gov", "ls", dirName], capture_output=True)
        listOfDirs = output.stdout.decode("utf-8").split("\n")
        for eachDir in listOfDirs:
            if eachDir != '':
                fullPath = dirName + eachDir
                output = subprocess.run(["eos", "root://cmseos.fnal.gov", "ls", fullPath], capture_output=True)
                listOfFiles = output.stdout.decode("utf-8").split("\n")
                for fileName in listOfFiles:
                    if fileName != '':
                        f.write("root://cmseos.fnal.gov/" + fullPath + '/' + fileName + '\n')
    else:
        print("Directory not found... Failed to write.")