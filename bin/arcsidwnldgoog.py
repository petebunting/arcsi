#! /usr/bin/env python

"""
Module that contains the ARSCI command to download imagery using gsutils
"""

############################################################################
#  arcsidwnldgoog.py
#
#  Copyright 2017 ARCSI.
#
#  ARCSI: 'Atmospheric and Radiometric Correction of Satellite Imagery'
#
#  ARCSI is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  ARCSI is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with ARCSI.  If not, see <http://www.gnu.org/licenses/>.
#
#
# Purpose:  A script to download imagery from Google given list of files.
#
# Author: Pete Bunting
# Email: pfb@aber.ac.uk
# Date: 14/07/2017
# Version: 1.0
#
# History:
# Version 1.0 - Created.
#
############################################################################

# Import updated print function into python 2.7
from __future__ import print_function
# Import updated division operator into python 2.7
from __future__ import division
# Import the python os.path module
import os.path
# Import the python Argument parser
import argparse
# Import python subprocess module
import subprocess
# Import the sys module
import sys

def readTextFile2List(file):
    """
    Read a text file into a list where each line 
    is an element in the list.
    """
    outList = []
    try:
        dataFile = open(file, 'r')
        for line in dataFile:
            line = line.strip()
            if line != "":
                outList.append(line)
        dataFile.close()
    except Exception as e:
        raise e
    return outList

def writeList2File(dataList, outFile):
    """
    Write a list a text file, one line per item.
    """
    try:
        f = open(outFile, 'w')
        for item in dataList:
           f.write(str(item)+'\n')
        f.flush()
        f.close()
    except Exception as e:
        raise e

def runGoogleImgDwbld(inFileDwnLst, outDIR, outFailLst=None, multiDwn=False, overwrite=False):

    fileLst = readTextFile2List(inFileDwnLst)
    failsLst = []

    multiStr = ''
    if multiDwn:
        multiStr = '-m'

    for file in fileLst:
        outFileName = os.path.basename(file)
        print("Processing " + outFileName)
        if overwrite or (not os.path.exists(os.path.join(outDIR, outFileName))):
            cmd = "gsutil "+multiStr+" cp -r " + file + " " + outDIR
            try:
                subprocess.call(cmd, shell=True)
            except OSError as e:
                failsLst.append(file)
                print("Error: {}".format(e), file=sys.stderr)
            except Exception as e:
                failsLst.append(file)
                print("Error: {}".format(e), file=sys.stderr)
        else:
            print("\tAlready Downloaded...")

        if not outFailLst is None:
            writeList2File(failsLst, outFailLst)


if __name__ == '__main__':
    """
    The command line user interface to ARCSI to download a list of files from google.
    """
    parser = argparse.ArgumentParser(prog='arcsidwnldgoog.py',
                                    description='''ARSCI command to download imagery from Google bucket''',
                                    epilog='''A tool to download imagery from a Google Bucket.''')

    parser.add_argument("-i", "--input", type=str, required=True, help='''Input file which lists gs:// paths to be downloaded.''')
    parser.add_argument("-o", "--outpath", type=str, required=True, help='''Output directory path where downloads to be downloaded to on your system.''')
    parser.add_argument("--fails", type=str, help='''Output file which lists any downloads which fail.''')
    parser.add_argument("--multi", action='store_true', default=False, help='''Adds -m option to the gsutil download command.''')
    parser.add_argument("--overwrite", action='store_true', default=False, help='''Redownloads and overwrites existing images, otherwise files which exist are not redownloaded.''')

    # Call the parser to parse the arguments.
    args = parser.parse_args()

    runGoogleImgDwbld(args.input, args.outpath, args.fails, args.multi, args.overwrite)

