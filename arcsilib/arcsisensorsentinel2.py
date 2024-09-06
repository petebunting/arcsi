"""
Module that contains the ARCSISentinel2Sensor class.
"""
############################################################################
#  arcsisensorsentinel2.py
#
#  Copyright 2016 ARCSI.
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
# Purpose:  A class for read the Sentinel-2 sensor header file and applying
#           the pre-processing operations within ARCSI to the Sentinel-2
#           datasets.
#
# Author: Pete Bunting
# Email: pfb@aber.ac.uk
# Date: 09/07/2016
# Version: 1.0
#
# History:
# Version 1.0 - Created.
#
############################################################################

import collections
import datetime
import glob
import math
import os
import shutil
import statistics
import xml.etree.ElementTree as ET

import fmask.config
import fmask.fmask
import numpy
import Py6S
import rsgislib.imagecalc
import rsgislib.imagecalibration
import rsgislib.imagemorphology
import rsgislib.imageutils
import rsgislib.rastergis
import rsgislib.segmentation
import rsgislib.tools.geometrytools
import rsgislib.tools.utils
from osgeo import gdal, osr
from rios import rat

import arcsilib.arcsiutils

from .arcsiexception import ARCSIException
from .arcsisensor import ARCSIAbstractSensor


class ARCSISen2SpectralBandObj(object):
    """
    This is a class to store the information associated with a Sentinel-2 image band.
    """

    def __init__(
        self,
        phyBandName=None,
        bandID=None,
        imgRes=0.0,
        wvLenMin=0.0,
        wvLenMax=0.0,
        wvLenCen=0.0,
        respFuncStep=1.0,
        respFunc=None,
    ):
        if respFunc is None:
            respFunc = list()
        self.phyBandName = phyBandName
        self.bandID = bandID
        self.imgRes = imgRes
        self.wvLenMin = wvLenMin
        self.wvLenMax = wvLenMax
        self.wvLenCen = wvLenCen
        self.respFuncStep = respFuncStep
        self.respFunc = respFunc
        self.respFunc6S = None
        self.wvLenMin6S = 0.0
        self.wvLenMax6S = 0.0


class ARCSISentinel2Sensor(ARCSIAbstractSensor):
    """
    A class which represents the Sentinel-2 sensor to read
    header parameters and apply data processing operations.
    """

    def __init__(self, debugMode, inputImage):
        ARCSIAbstractSensor.__init__(self, debugMode, inputImage)
        self.sensor = "SEN2"
        self.sen2FileBaseDIR = ""
        self.processingLevel = ""
        self.processingBaseline = ""
        self.spacecraftName = ""
        self.dataTakeType = ""
        self.orbitNumber = ""
        self.orbitDirection = ""

        self.sen2ImgB01 = ""
        self.sen2ImgB02 = ""
        self.sen2ImgB03 = ""
        self.sen2ImgB04 = ""
        self.sen2ImgB05 = ""
        self.sen2ImgB06 = ""
        self.sen2ImgB07 = ""
        self.sen2ImgB8A = ""
        self.sen2ImgB08 = ""
        self.sen2ImgB09 = ""
        self.sen2ImgB10 = ""
        self.sen2ImgB11 = ""
        self.sen2ImgB12 = ""
        self.sen2ImgTCI = ""

        # File paths/names for resampled images
        self.sen2ImgB02_20m = None
        self.sen2ImgB03_20m = None
        self.sen2ImgB04_20m = None
        self.sen2ImgB08_20m = None
        self.sen2ImgB05_10m = None
        self.sen2ImgB06_10m = None
        self.sen2ImgB07_10m = None
        self.sen2ImgB8A_10m = None
        self.sen2ImgB11_10m = None
        self.sen2ImgB12_10m = None

        self.inNoDataVal = 0
        self.inSatDataVal = 65535
        self.quantificationVal = 10000

        self.physicalGain_B0 = 0.0
        self.physicalGain_B1 = 0.0
        self.physicalGain_B2 = 0.0
        self.physicalGain_B3 = 0.0
        self.physicalGain_B4 = 0.0
        self.physicalGain_B5 = 0.0
        self.physicalGain_B6 = 0.0
        self.physicalGain_B7 = 0.0
        self.physicalGain_B8 = 0.0
        self.physicalGain_B9 = 0.0
        self.physicalGain_B10 = 0.0
        self.physicalGain_B11 = 0.0
        self.physicalGain_B12 = 0.0

        self.earthSunDist_U = 0.0
        self.esun_B0 = 0.0
        self.esun_B1 = 0.0
        self.esun_B2 = 0.0
        self.esun_B3 = 0.0
        self.esun_B4 = 0.0
        self.esun_B5 = 0.0
        self.esun_B6 = 0.0
        self.esun_B7 = 0.0
        self.esun_B8 = 0.0
        self.esun_B9 = 0.0
        self.esun_B10 = 0.0
        self.esun_B11 = 0.0
        self.esun_B12 = 0.0

        self.ratiometric_offs_B0 = 0.0
        self.ratiometric_offs_B1 = 0.0
        self.ratiometric_offs_B2 = 0.0
        self.ratiometric_offs_B3 = 0.0
        self.ratiometric_offs_B4 = 0.0
        self.ratiometric_offs_B5 = 0.0
        self.ratiometric_offs_B6 = 0.0
        self.ratiometric_offs_B7 = 0.0
        self.ratiometric_offs_B8 = 0.0
        self.ratiometric_offs_B9 = 0.0
        self.ratiometric_offs_B10 = 0.0
        self.ratiometric_offs_B11 = 0.0
        self.ratiometric_offs_B12 = 0.0

        self.specBandInfo = dict()

        self.uniqueTileID = ""

        self.resampleTo20m = False

        self.epsgCode = 0

        self.imgIntScaleFactor = 1000

    ################################
    #### CODE FROM PYTHON FMASK ####
    ################################
    @staticmethod
    def makeValueArray(valuesListNode):
        """
        Take a <Values_List> node from the XML, and return an array of the values contained
        within it. This will be a 2-d numpy array of float32 values (should I pass the dtype in??)

        """
        valuesList = valuesListNode.findall("VALUES")
        vals = []
        for valNode in valuesList:
            text = valNode.text
            vals.append([numpy.float32(x) for x in text.strip().split()])
        return numpy.array(vals)

    def buildViewAngleArr(self, viewingAngleNodeList, angleName):
        """
        Build up the named viewing angle array from the various detector strips given as
        separate arrays. I don't really understand this, and may need to re-write it once
        I have worked it out......

        The angleName is one of 'Zenith' or 'Azimuth'.
        Returns a dictionary of 2-d arrays, keyed by the bandId string.
        """
        angleArrDict = {}
        for viewingAngleNode in viewingAngleNodeList:
            bandId = viewingAngleNode.attrib["bandId"]
            angleNode = viewingAngleNode.find(angleName)
            angleArr = self.makeValueArray(angleNode.find("Values_List"))
            if bandId not in angleArrDict:
                angleArrDict[bandId] = angleArr
            else:
                mask = ~numpy.isnan(angleArr)
                angleArrDict[bandId][mask] = angleArr[mask]
        return angleArrDict

    ################################

    def extractHeaderParameters(self, inputHeader, wktStr):
        """
        Understands and parses the Sentinel-2 xml header file
        """
        try:

            inputHeader = os.path.abspath(inputHeader)
            self.headerFileName = os.path.split(inputHeader)[1]
            self.sen2FileBaseDIR = os.path.split(inputHeader)[0]

            tree = ET.parse(inputHeader)
            root = tree.getroot()

            try:
                schemaLocation = root.attrib.get('{http://www.w3.org/2001/XMLSchema-instance}schemaLocation')
                psd_version = int(schemaLocation.split(".")[0].split("-")[-1])
                tag_search = f"{{https://psd-{psd_version}.sentinel2.eo.esa.int/PSD/User_Product_Level-1C.xsd}}General_Info"
                generalInfoTag = root.find(tag_search)

                if generalInfoTag is None:
                    raise ARCSIException(f'Unable to find General_Info tag using search: {tag_search}')

                if psd_version == 13:
                    hdrFileVersion = "psd13"
                elif psd_version >= 14:
                    # force version15 to version14 (version15 is introduced in baseline 5.11 on 2024-07-23)
                    hdrFileVersion = "psd14"
                else:
                    raise ARCSIException(
                    f"PSD version = {psd_version}. Expected versions 13 or above"
                )

            except Exception as error:
                raise ARCSIException(
                    f"Problem parsing xml 'General_Info' tag. Message = {error}"
                )

            productInfoTag = generalInfoTag.find("Product_Info")
            if productInfoTag == None:
                raise ARCSIException(
                    "Cannot open 'General_Info' section 'Product_Info' - is this really a Sentinel-2 image file?"
                )

            productType = productInfoTag.find("PRODUCT_TYPE").text.strip()
            if (productType == None) or (productType == ""):
                raise ARCSIException(
                    "Cannot find the product type - is this really a Sentinel-2 image file?"
                )
            elif productType != "S2MSI1C":
                raise ARCSIException("Expecting Sentinel-2 product type 'S2MSI1C'")

            self.processingLevel = productInfoTag.find("PROCESSING_LEVEL").text.strip()
            self.processingBaseline = productInfoTag.find(
                "PROCESSING_BASELINE"
            ).text.strip()

            genTimeStr = productInfoTag.find("GENERATION_TIME").text.strip()
            genTimeStr = genTimeStr.replace("Z", "")
            self.generationTime = datetime.datetime.strptime(
                genTimeStr, "%Y-%m-%dT%H:%M:%S.%f"
            )

            datatakeTag = productInfoTag.find("Datatake")
            if datatakeTag == None:
                raise ARCSIException("Cannot open 'Product_Info' section 'Datatake'")

            acqTimeStr = datatakeTag.find("DATATAKE_SENSING_START").text.strip()
            acqTimeStr = acqTimeStr.replace("Z", "")
            self.acquisitionTime = datetime.datetime.strptime(
                acqTimeStr, "%Y-%m-%dT%H:%M:%S.%f"
            )

            self.spacecraftName = datatakeTag.find("SPACECRAFT_NAME").text.strip()
            self.dataTakeType = datatakeTag.find("DATATAKE_TYPE").text.strip()
            self.orbitNumber = datatakeTag.find("SENSING_ORBIT_NUMBER").text.strip()
            self.orbitDirection = datatakeTag.find(
                "SENSING_ORBIT_DIRECTION"
            ).text.strip()

            haveUniqueTileID = False
            mixedPSD14PSD13 = False
            if hdrFileVersion == "psd14":
                prodURI = productInfoTag.find("PRODUCT_URI").text.strip()
                prodURIList = prodURI.split("_")
                if (
                    ("OPER" in prodURIList)
                    and ("REQ" in prodURIList)
                    and ("PRDDWN" in prodURIList)
                ):
                    haveUniqueTileID = False
                    mixedPSD14PSD13 = True
                elif len(prodURIList) >= 6:
                    self.uniqueTileID = prodURIList[5]
                    haveUniqueTileID = True
                else:
                    haveUniqueTileID = False
                    mixedPSD14PSD13 = True

            # Get the input granule information (e.g., bands).
            granulesTagsLst = list()
            if hdrFileVersion == "psd14":
                granuleListTag = productInfoTag.find("Product_Organisation").find(
                    "Granule_List"
                )
                for granuleLstChild in granuleListTag:
                    if granuleLstChild.tag == "Granule":
                        granulesTagsLst.append(granuleLstChild)
            elif hdrFileVersion == "psd13":
                granuleListTags = productInfoTag.find("Product_Organisation")
                for granulesLstTag in granuleListTags:
                    if granulesLstTag.tag == "Granule_List":
                        granulesTagsLst.append(granulesLstTag.find("Granules"))
            else:
                raise ARCSIException(
                    "Do not recognised header format: '" + hdrFileVersion + "'"
                )

            granuleIdentifier = ""
            if len(granulesTagsLst) == 0:
                raise ARCSIException(
                    "Could not find any granules within the header file. (Header: "
                    + hdrFileVersion
                    + ")"
                )
            elif (hdrFileVersion == "psd13") or (
                (len(granulesTagsLst) != 1)
                and (hdrFileVersion == "psd14")
                and mixedPSD14PSD13
            ):
                granulesDIR = os.path.join(self.sen2FileBaseDIR, "GRANULE")
                files = os.listdir(path=granulesDIR)
                granuleFileTags = list()
                for file in files:
                    if os.path.isdir(os.path.join(granulesDIR, file)):
                        for granulesTmpTag in granulesTagsLst:
                            if (
                                granulesTmpTag.attrib["granuleIdentifier"].strip()
                                == file
                            ):
                                granuleFileTags.append(granulesTmpTag)
                if len(granuleFileTags) == 1:
                    granuleTag = granuleFileTags[0]
                    granuleIdentifier = granuleTag.attrib["granuleIdentifier"].strip()
                    self.uniqueTileID = granuleIdentifier.split("_")[9]
                    haveUniqueTileID = True
                else:
                    raise ARCSIException(
                        "Only expecting a single granule within the file... The input image you have provided is not supported by ARCSI - please report so we can add support."
                    )
            elif (len(granulesTagsLst) != 1) and (hdrFileVersion == "psd14"):
                raise ARCSIException(
                    "Only expecting a single granule within the file... The input image you have provided is not supported by ARCSI - please report so we can add support."
                )
            else:
                granuleTag = granulesTagsLst[0]

            if (hdrFileVersion == "psd14") and (not mixedPSD14PSD13):
                for granuleChild in granuleTag:
                    if granuleChild.tag == "IMAGE_FILE":
                        imgFile = granuleChild.text.strip()
                        tmpFiles = glob.glob(
                            os.path.join(self.sen2FileBaseDIR, imgFile + "*.jp2")
                        )
                        if len(tmpFiles) == 1:
                            if "B01" in imgFile:
                                self.sen2ImgB01 = tmpFiles[0]
                            elif "B02" in imgFile:
                                self.sen2ImgB02 = tmpFiles[0]
                            elif "B03" in imgFile:
                                self.sen2ImgB03 = tmpFiles[0]
                            elif "B04" in imgFile:
                                self.sen2ImgB04 = tmpFiles[0]
                            elif "B05" in imgFile:
                                self.sen2ImgB05 = tmpFiles[0]
                            elif "B06" in imgFile:
                                self.sen2ImgB06 = tmpFiles[0]
                            elif "B07" in imgFile:
                                self.sen2ImgB07 = tmpFiles[0]
                            elif "B8A" in imgFile:
                                self.sen2ImgB8A = tmpFiles[0]
                            elif "B08" in imgFile:
                                self.sen2ImgB08 = tmpFiles[0]
                            elif "B09" in imgFile:
                                self.sen2ImgB09 = tmpFiles[0]
                            elif "B10" in imgFile:
                                self.sen2ImgB10 = tmpFiles[0]
                            elif "B11" in imgFile:
                                self.sen2ImgB11 = tmpFiles[0]
                            elif "B12" in imgFile:
                                self.sen2ImgB12 = tmpFiles[0]
                            elif "TCI" in imgFile:
                                self.sen2ImgTCI = tmpFiles[0]
                            else:
                                raise ARCSIException(
                                    "Could not associated image file with an expected image band: "
                                    + imgFile
                                )
                        else:
                            raise ARCSIException(
                                "Could not file image file for: " + imgFile
                            )
                    else:
                        raise ARCSIException("There is no tag 'IMAGE_FILE'.")
            elif (hdrFileVersion == "psd13") or (
                (hdrFileVersion == "psd14") and mixedPSD14PSD13
            ):
                imgsPath = os.path.join(
                    os.path.join(
                        os.path.join(self.sen2FileBaseDIR, "GRANULE"), granuleIdentifier
                    ),
                    "IMG_DATA",
                )
                for granuleChild in granuleTag:
                    if granuleChild.tag == "IMAGE_ID":
                        imgFile = granuleChild.text.strip()
                        tmpFiles = glob.glob(os.path.join(imgsPath, imgFile + "*.jp2"))
                        if len(tmpFiles) == 1:
                            if "B01" in imgFile:
                                self.sen2ImgB01 = tmpFiles[0]
                            elif "B02" in imgFile:
                                self.sen2ImgB02 = tmpFiles[0]
                            elif "B03" in imgFile:
                                self.sen2ImgB03 = tmpFiles[0]
                            elif "B04" in imgFile:
                                self.sen2ImgB04 = tmpFiles[0]
                            elif "B05" in imgFile:
                                self.sen2ImgB05 = tmpFiles[0]
                            elif "B06" in imgFile:
                                self.sen2ImgB06 = tmpFiles[0]
                            elif "B07" in imgFile:
                                self.sen2ImgB07 = tmpFiles[0]
                            elif "B8A" in imgFile:
                                self.sen2ImgB8A = tmpFiles[0]
                            elif "B08" in imgFile:
                                self.sen2ImgB08 = tmpFiles[0]
                            elif "B09" in imgFile:
                                self.sen2ImgB09 = tmpFiles[0]
                            elif "B10" in imgFile:
                                self.sen2ImgB10 = tmpFiles[0]
                            elif "B11" in imgFile:
                                self.sen2ImgB11 = tmpFiles[0]
                            elif "B12" in imgFile:
                                self.sen2ImgB12 = tmpFiles[0]
                            elif "TCI" in imgFile:
                                self.sen2ImgTCI = tmpFiles[0]
                            else:
                                raise ARCSIException(
                                    "Could not associated image file with an expected image band: "
                                    + imgFile
                                )
                        else:
                            raise ARCSIException(
                                "Could not file image file for: " + imgFile
                            )
            else:
                raise ARCSIException(
                    "Do not recognised header format: '" + hdrFileVersion + "'"
                )

            productImgCharTag = generalInfoTag.find("Product_Image_Characteristics")
            if productImgCharTag == None:
                raise ARCSIException(
                    "Cannot open 'General_Info' section 'Product_Image_Characteristics'"
                )

            specialValsTagsLst = list()
            quantificationValTag = None
            reflectanceConversionTag = None
            spectralInfoListTag = None
            radiometricOffsetListTag = None
            physicalGainsTagsLst = list()
            for prodImgCharChild in productImgCharTag:
                if prodImgCharChild.tag == "Special_Values":
                    specialValsTagsLst.append(prodImgCharChild)
                elif prodImgCharChild.tag == "QUANTIFICATION_VALUE":
                    quantificationValTag = prodImgCharChild
                elif prodImgCharChild.tag == "Reflectance_Conversion":
                    reflectanceConversionTag = prodImgCharChild
                elif prodImgCharChild.tag == "Spectral_Information_List":
                    spectralInfoListTag = prodImgCharChild
                elif prodImgCharChild.tag == "PHYSICAL_GAINS":
                    physicalGainsTagsLst.append(prodImgCharChild)
                elif prodImgCharChild.tag == "Radiometric_Offset_List":
                    radiometricOffsetListTag = prodImgCharChild

            # Parse out the input no data and saturation value.
            for tag in specialValsTagsLst:
                if tag.find("SPECIAL_VALUE_TEXT").text.strip() == "NODATA":
                    self.inNoDataVal = rsgislib.tools.utils.str_to_int(
                        tag.find("SPECIAL_VALUE_INDEX").text.strip()
                    )
                elif tag.find("SPECIAL_VALUE_TEXT").text.strip() == "SATURATED":
                    self.inSatDataVal = rsgislib.tools.utils.str_to_int(
                        tag.find("SPECIAL_VALUE_INDEX").text.strip()
                    )

            # Get the Quantification value.
            self.quantificationVal = rsgislib.tools.utils.str_to_float(
                quantificationValTag.text.strip()
            )

            # Get radiometric offsets (is provided; introduced 79.90)
            if radiometricOffsetListTag is not None:
                for radiometricOffsetTag in radiometricOffsetListTag:
                    if radiometricOffsetTag.attrib["band_id"] == "0":
                        self.ratiometric_offs_B0 = rsgislib.tools.utils.str_to_float(
                            radiometricOffsetTag.text.strip()
                        )
                    elif radiometricOffsetTag.attrib["band_id"] == "1":
                        self.ratiometric_offs_B1 = rsgislib.tools.utils.str_to_float(
                            radiometricOffsetTag.text.strip()
                        )
                    elif radiometricOffsetTag.attrib["band_id"] == "2":
                        self.ratiometric_offs_B2 = rsgislib.tools.utils.str_to_float(
                            radiometricOffsetTag.text.strip()
                        )
                    elif radiometricOffsetTag.attrib["band_id"] == "3":
                        self.ratiometric_offs_B3 = rsgislib.tools.utils.str_to_float(
                            radiometricOffsetTag.text.strip()
                        )
                    elif radiometricOffsetTag.attrib["band_id"] == "4":
                        self.ratiometric_offs_B4 = rsgislib.tools.utils.str_to_float(
                            radiometricOffsetTag.text.strip()
                        )
                    elif radiometricOffsetTag.attrib["band_id"] == "5":
                        self.ratiometric_offs_B5 = rsgislib.tools.utils.str_to_float(
                            radiometricOffsetTag.text.strip()
                        )
                    elif radiometricOffsetTag.attrib["band_id"] == "6":
                        self.ratiometric_offs_B6 = rsgislib.tools.utils.str_to_float(
                            radiometricOffsetTag.text.strip()
                        )
                    elif radiometricOffsetTag.attrib["band_id"] == "7":
                        self.ratiometric_offs_B7 = rsgislib.tools.utils.str_to_float(
                            radiometricOffsetTag.text.strip()
                        )
                    elif radiometricOffsetTag.attrib["band_id"] == "8":
                        self.ratiometric_offs_B8 = rsgislib.tools.utils.str_to_float(
                            radiometricOffsetTag.text.strip()
                        )
                    elif radiometricOffsetTag.attrib["band_id"] == "9":
                        self.ratiometric_offs_B9 = rsgislib.tools.utils.str_to_float(
                            radiometricOffsetTag.text.strip()
                        )
                    elif radiometricOffsetTag.attrib["band_id"] == "10":
                        self.ratiometric_offs_B10 = rsgislib.tools.utils.str_to_float(
                            radiometricOffsetTag.text.strip()
                        )
                    elif radiometricOffsetTag.attrib["band_id"] == "11":
                        self.ratiometric_offs_B11 = rsgislib.tools.utils.str_to_float(
                            radiometricOffsetTag.text.strip()
                        )
                    elif radiometricOffsetTag.attrib["band_id"] == "12":
                        self.ratiometric_offs_B12 = rsgislib.tools.utils.str_to_float(
                            radiometricOffsetTag.text.strip()
                        )

            # Get the Physical Gain values.
            for physGainTag in physicalGainsTagsLst:
                if physGainTag.attrib["bandId"] == "0":
                    self.physicalGain_B0 = rsgislib.tools.utils.str_to_float(
                        physGainTag.text.strip()
                    )
                elif physGainTag.attrib["bandId"] == "1":
                    self.physicalGain_B1 = rsgislib.tools.utils.str_to_float(
                        physGainTag.text.strip()
                    )
                elif physGainTag.attrib["bandId"] == "2":
                    self.physicalGain_B2 = rsgislib.tools.utils.str_to_float(
                        physGainTag.text.strip()
                    )
                elif physGainTag.attrib["bandId"] == "3":
                    self.physicalGain_B3 = rsgislib.tools.utils.str_to_float(
                        physGainTag.text.strip()
                    )
                elif physGainTag.attrib["bandId"] == "4":
                    self.physicalGain_B4 = rsgislib.tools.utils.str_to_float(
                        physGainTag.text.strip()
                    )
                elif physGainTag.attrib["bandId"] == "5":
                    self.physicalGain_B5 = rsgislib.tools.utils.str_to_float(
                        physGainTag.text.strip()
                    )
                elif physGainTag.attrib["bandId"] == "6":
                    self.physicalGain_B6 = rsgislib.tools.utils.str_to_float(
                        physGainTag.text.strip()
                    )
                elif physGainTag.attrib["bandId"] == "7":
                    self.physicalGain_B7 = rsgislib.tools.utils.str_to_float(
                        physGainTag.text.strip()
                    )
                elif physGainTag.attrib["bandId"] == "8":
                    self.physicalGain_B8 = rsgislib.tools.utils.str_to_float(
                        physGainTag.text.strip()
                    )
                elif physGainTag.attrib["bandId"] == "9":
                    self.physicalGain_B9 = rsgislib.tools.utils.str_to_float(
                        physGainTag.text.strip()
                    )
                elif physGainTag.attrib["bandId"] == "10":
                    self.physicalGain_B10 = rsgislib.tools.utils.str_to_float(
                        physGainTag.text.strip()
                    )
                elif physGainTag.attrib["bandId"] == "11":
                    self.physicalGain_B11 = rsgislib.tools.utils.str_to_float(
                        physGainTag.text.strip()
                    )
                elif physGainTag.attrib["bandId"] == "12":
                    self.physicalGain_B12 = rsgislib.tools.utils.str_to_float(
                        physGainTag.text.strip()
                    )

            # Get Reflectance Conversion info
            self.earthSunDist_U = rsgislib.tools.utils.str_to_float(
                reflectanceConversionTag.find("U").text.strip()
            )
            for solarIrrTag in reflectanceConversionTag.find("Solar_Irradiance_List"):
                if solarIrrTag.attrib["bandId"] == "0":
                    self.esun_B0 = rsgislib.tools.utils.str_to_float(
                        solarIrrTag.text.strip()
                    )
                elif solarIrrTag.attrib["bandId"] == "1":
                    self.esun_B1 = rsgislib.tools.utils.str_to_float(
                        solarIrrTag.text.strip()
                    )
                elif solarIrrTag.attrib["bandId"] == "2":
                    self.esun_B2 = rsgislib.tools.utils.str_to_float(
                        solarIrrTag.text.strip()
                    )
                elif solarIrrTag.attrib["bandId"] == "3":
                    self.esun_B3 = rsgislib.tools.utils.str_to_float(
                        solarIrrTag.text.strip()
                    )
                elif solarIrrTag.attrib["bandId"] == "4":
                    self.esun_B4 = rsgislib.tools.utils.str_to_float(
                        solarIrrTag.text.strip()
                    )
                elif solarIrrTag.attrib["bandId"] == "5":
                    self.esun_B5 = rsgislib.tools.utils.str_to_float(
                        solarIrrTag.text.strip()
                    )
                elif solarIrrTag.attrib["bandId"] == "6":
                    self.esun_B6 = rsgislib.tools.utils.str_to_float(
                        solarIrrTag.text.strip()
                    )
                elif solarIrrTag.attrib["bandId"] == "7":
                    self.esun_B7 = rsgislib.tools.utils.str_to_float(
                        solarIrrTag.text.strip()
                    )
                elif solarIrrTag.attrib["bandId"] == "8":
                    self.esun_B8 = rsgislib.tools.utils.str_to_float(
                        solarIrrTag.text.strip()
                    )
                elif solarIrrTag.attrib["bandId"] == "9":
                    self.esun_B9 = rsgislib.tools.utils.str_to_float(
                        solarIrrTag.text.strip()
                    )
                elif solarIrrTag.attrib["bandId"] == "10":
                    self.esun_B10 = rsgislib.tools.utils.str_to_float(
                        solarIrrTag.text.strip()
                    )
                elif solarIrrTag.attrib["bandId"] == "11":
                    self.esun_B11 = rsgislib.tools.utils.str_to_float(
                        solarIrrTag.text.strip()
                    )
                elif solarIrrTag.attrib["bandId"] == "12":
                    self.esun_B12 = rsgislib.tools.utils.str_to_float(
                        solarIrrTag.text.strip()
                    )

            if (hdrFileVersion == "psd14") and (spectralInfoListTag is not None):
                # Get Spectral Info
                for specInfoTag in spectralInfoListTag:
                    phyBandName = specInfoTag.attrib["physicalBand"]
                    specBandObj = ARCSISen2SpectralBandObj()
                    specBandObj.phyBandName = phyBandName
                    specBandObj.bandID = specInfoTag.attrib["bandId"]
                    specBandObj.imgRes = rsgislib.tools.utils.str_to_float(
                        specInfoTag.find("RESOLUTION").text.strip()
                    )
                    specBandObj.wvLenMin = rsgislib.tools.utils.str_to_float(
                        specInfoTag.find("Wavelength").find("MIN").text.strip()
                    )
                    specBandObj.wvLenMax = rsgislib.tools.utils.str_to_float(
                        specInfoTag.find("Wavelength").find("MAX").text.strip()
                    )
                    specBandObj.wvLenCen = rsgislib.tools.utils.str_to_float(
                        specInfoTag.find("Wavelength").find("CENTRAL").text.strip()
                    )
                    specBandObj.respFuncStep = rsgislib.tools.utils.str_to_float(
                        specInfoTag.find("Spectral_Response").find("STEP").text.strip()
                    )
                    specResStrVals = (
                        specInfoTag.find("Spectral_Response")
                        .find("VALUES")
                        .text.strip()
                        .split(" ")
                    )
                    specResVals = list()
                    for strVal in specResStrVals:
                        specResVals.append(rsgislib.tools.utils.str_to_float(strVal))
                    specBandObj.respFunc = specResVals
                    self.specBandInfo[phyBandName] = specBandObj
            else:
                # B01
                specBandObj = ARCSISen2SpectralBandObj()
                specBandObj.phyBandName = "B1"
                specBandObj.bandID = "0"
                specBandObj.imgRes = 60
                specBandObj.wvLenMin = 430
                specBandObj.wvLenMax = 457
                specBandObj.wvLenCen = 443.9
                specBandObj.respFuncStep = 1
                specBandObj.respFunc = [
                    0.01522444,
                    0.06669758,
                    0.19425897,
                    0.35395736,
                    0.45648857,
                    0.50759455,
                    0.54750739,
                    0.58419244,
                    0.61012868,
                    0.64603585,
                    0.69458246,
                    0.74037505,
                    0.78703023,
                    0.85862712,
                    0.94458791,
                    0.9928916,
                    1,
                    0.99055275,
                    0.97282606,
                    0.95596914,
                    0.95429069,
                    0.91888272,
                    0.72055356,
                    0.38639386,
                    0.14531035,
                    0.05161255,
                    0.01738704,
                    0.00029585,
                ]
                self.specBandInfo["B1"] = specBandObj

                # B02
                specBandObj = ARCSISen2SpectralBandObj()
                specBandObj.phyBandName = "B2"
                specBandObj.bandID = "1"
                specBandObj.imgRes = 10
                specBandObj.wvLenMin = 440
                specBandObj.wvLenMax = 538
                specBandObj.wvLenCen = 496.6
                specBandObj.respFuncStep = 1
                specBandObj.respFunc = [
                    0.00119988,
                    0.00201397,
                    0.00258793,
                    0.00271734,
                    0.00271858,
                    0.002053,
                    0.00324912,
                    0.0021993,
                    0.00277292,
                    0.00311194,
                    0.00234723,
                    0.00374245,
                    0.0028408,
                    0.00304821,
                    0.00604983,
                    0.00894596,
                    0.01953246,
                    0.03875845,
                    0.08374839,
                    0.17565347,
                    0.29129289,
                    0.36347223,
                    0.3811347,
                    0.38419864,
                    0.39176673,
                    0.39862405,
                    0.40894049,
                    0.42354641,
                    0.4485657,
                    0.4811418,
                    0.50498541,
                    0.52293008,
                    0.52892822,
                    0.53366,
                    0.53242234,
                    0.53311303,
                    0.53655971,
                    0.54232711,
                    0.55667534,
                    0.57791322,
                    0.60145975,
                    0.6156357,
                    0.62060573,
                    0.61270938,
                    0.59482968,
                    0.57420278,
                    0.55609253,
                    0.5440646,
                    0.54004284,
                    0.5517318,
                    0.56998769,
                    0.59684728,
                    0.63205242,
                    0.67244298,
                    0.71093613,
                    0.73748447,
                    0.75709994,
                    0.76697185,
                    0.77176039,
                    0.77883444,
                    0.78683055,
                    0.79421954,
                    0.80824012,
                    0.82348832,
                    0.83743831,
                    0.84485726,
                    0.84716089,
                    0.83974417,
                    0.82502148,
                    0.8036499,
                    0.78544282,
                    0.76973497,
                    0.7598602,
                    0.76337273,
                    0.77981251,
                    0.80847605,
                    0.84947272,
                    0.90112566,
                    0.95456662,
                    0.98736039,
                    1,
                    0.98609155,
                    0.90770989,
                    0.72315884,
                    0.47814326,
                    0.28641509,
                    0.16955089,
                    0.10257285,
                    0.06498784,
                    0.04106167,
                    0.02503855,
                    0.01307564,
                    0.00257814,
                    0.00108051,
                    0.00030609,
                    0.00043924,
                    0.00044121,
                ]
                self.specBandInfo["B2"] = specBandObj

                # B03
                specBandObj = ARCSISen2SpectralBandObj()
                specBandObj.phyBandName = "B3"
                specBandObj.bandID = "2"
                specBandObj.imgRes = 10
                specBandObj.wvLenMin = 537
                specBandObj.wvLenMax = 582
                specBandObj.wvLenCen = 560
                specBandObj.respFuncStep = 1
                specBandObj.respFunc = [
                    0.00080152,
                    0.01631966,
                    0.03749604,
                    0.08021834,
                    0.16857673,
                    0.33961135,
                    0.57045802,
                    0.74395167,
                    0.8255379,
                    0.86623109,
                    0.88713486,
                    0.89063153,
                    0.87743881,
                    0.85952176,
                    0.84272738,
                    0.83271245,
                    0.83091319,
                    0.8429466,
                    0.86557037,
                    0.89523547,
                    0.93204973,
                    0.96550034,
                    0.99001699,
                    1,
                    0.99850933,
                    0.98241577,
                    0.94879561,
                    0.90893224,
                    0.87016848,
                    0.83868631,
                    0.8133992,
                    0.79225145,
                    0.7842798,
                    0.78830002,
                    0.80532973,
                    0.82861237,
                    0.84453213,
                    0.85667749,
                    0.85654311,
                    0.79885992,
                    0.62453426,
                    0.38688244,
                    0.20018537,
                    0.09831467,
                    0.04284073,
                    0.01651146,
                ]
                self.specBandInfo["B3"] = specBandObj

                # B04
                specBandObj = ARCSISen2SpectralBandObj()
                specBandObj.phyBandName = "B4"
                specBandObj.bandID = 3
                specBandObj.imgRes = 10
                specBandObj.wvLenMin = 646
                specBandObj.wvLenMax = 684
                specBandObj.wvLenCen = 664.5
                specBandObj.respFuncStep = 1
                specBandObj.respFunc = [
                    0.00261427,
                    0.03462832,
                    0.15030251,
                    0.46548409,
                    0.81834707,
                    0.96554871,
                    0.98388489,
                    0.99687187,
                    1,
                    0.9955785,
                    0.99164257,
                    0.97772062,
                    0.93750282,
                    0.87465366,
                    0.81520176,
                    0.77787363,
                    0.7662682,
                    0.77666981,
                    0.80308737,
                    0.83262125,
                    0.8589057,
                    0.88527593,
                    0.91047688,
                    0.93604508,
                    0.95692399,
                    0.96878538,
                    0.9736139,
                    0.97172876,
                    0.96901499,
                    0.96568155,
                    0.96045441,
                    0.94488073,
                    0.88430524,
                    0.70624874,
                    0.42290429,
                    0.18976191,
                    0.06313289,
                    0.02061386,
                    0.0020257,
                ]
                self.specBandInfo["B4"] = specBandObj

                # B05
                specBandObj = ARCSISen2SpectralBandObj()
                specBandObj.phyBandName = "B5"
                specBandObj.bandID = "4"
                specBandObj.imgRes = 20
                specBandObj.wvLenMin = 694
                specBandObj.wvLenMax = 713
                specBandObj.wvLenCen = 703.9
                specBandObj.respFuncStep = 1
                specBandObj.respFunc = [
                    0.00118221,
                    0.04128719,
                    0.16781115,
                    0.47867704,
                    0.83408915,
                    0.98555238,
                    1,
                    0.99917704,
                    0.99301208,
                    0.98202139,
                    0.96500594,
                    0.94523647,
                    0.92390813,
                    0.90154471,
                    0.88461764,
                    0.86012379,
                    0.75605334,
                    0.52042972,
                    0.19640628,
                    0.03678278,
                ]
                self.specBandInfo["B5"] = specBandObj

                # B06
                specBandObj = ARCSISen2SpectralBandObj()
                specBandObj.phyBandName = "B6"
                specBandObj.bandID = "5"
                specBandObj.imgRes = 20
                specBandObj.wvLenMin = 731
                specBandObj.wvLenMax = 749
                specBandObj.wvLenCen = 740.2
                specBandObj.respFuncStep = 1
                specBandObj.respFunc = [
                    0.00528628,
                    0.08491265,
                    0.34549055,
                    0.75026111,
                    0.91998424,
                    0.91774468,
                    0.93414364,
                    0.95786657,
                    0.97589351,
                    0.98201154,
                    0.98159765,
                    0.99345282,
                    1,
                    0.98250656,
                    0.96245634,
                    0.85475636,
                    0.50661225,
                    0.13533181,
                    0.0134302,
                ]
                self.specBandInfo["B6"] = specBandObj

                # B07
                specBandObj = ARCSISen2SpectralBandObj()
                specBandObj.phyBandName = "B7"
                specBandObj.bandID = "6"
                specBandObj.imgRes = 20
                specBandObj.wvLenMin = 769
                specBandObj.wvLenMax = 797
                specBandObj.wvLenCen = 782.5
                specBandObj.respFuncStep = 1
                specBandObj.respFunc = [
                    0.00158775,
                    0.01471955,
                    0.06700855,
                    0.19944036,
                    0.42271848,
                    0.69391142,
                    0.89840316,
                    0.98314165,
                    0.99479749,
                    1,
                    0.99483279,
                    0.96447136,
                    0.90781386,
                    0.8464478,
                    0.80150314,
                    0.77808053,
                    0.77627582,
                    0.78832546,
                    0.79959911,
                    0.80136031,
                    0.79006668,
                    0.75603297,
                    0.67647373,
                    0.53577608,
                    0.36341065,
                    0.19325756,
                    0.07716074,
                    0.01971336,
                    0.00315275,
                ]
                self.specBandInfo["B7"] = specBandObj

                # B08
                specBandObj = ARCSISen2SpectralBandObj()
                specBandObj.phyBandName = "B8"
                specBandObj.bandID = "7"
                specBandObj.imgRes = 10
                specBandObj.wvLenMin = 760
                specBandObj.wvLenMax = 908
                specBandObj.wvLenCen = 835.1
                specBandObj.respFuncStep = 1
                specBandObj.respFunc = [
                    0.00135242,
                    0.00391616,
                    0.00044871,
                    0.00759275,
                    0.01905313,
                    0.03349108,
                    0.05649128,
                    0.0870686,
                    0.13235321,
                    0.20327639,
                    0.31387542,
                    0.44988941,
                    0.58726605,
                    0.71436889,
                    0.8181812,
                    0.90284514,
                    0.96067672,
                    0.99369744,
                    1,
                    0.98524291,
                    0.95844788,
                    0.93666123,
                    0.92594982,
                    0.93050611,
                    0.94139304,
                    0.95341007,
                    0.96218307,
                    0.9655653,
                    0.96296703,
                    0.95877093,
                    0.95087228,
                    0.94471788,
                    0.94260088,
                    0.94521458,
                    0.94468494,
                    0.94302291,
                    0.9363001,
                    0.92707231,
                    0.91511356,
                    0.90021968,
                    0.88081425,
                    0.86148256,
                    0.84257439,
                    0.82215879,
                    0.80140132,
                    0.7765823,
                    0.75539136,
                    0.73775889,
                    0.72215744,
                    0.70870534,
                    0.69854507,
                    0.6903735,
                    0.68251717,
                    0.68178973,
                    0.68302899,
                    0.67891416,
                    0.67639408,
                    0.67176564,
                    0.66600791,
                    0.66127505,
                    0.65915263,
                    0.65868891,
                    0.66436872,
                    0.67295613,
                    0.68563017,
                    0.7011901,
                    0.72062162,
                    0.74210801,
                    0.75925571,
                    0.77620597,
                    0.7835688,
                    0.78713055,
                    0.78702403,
                    0.7828085,
                    0.77539043,
                    0.7675732,
                    0.75848677,
                    0.74517599,
                    0.73227212,
                    0.71988842,
                    0.70601879,
                    0.69027923,
                    0.67538468,
                    0.66109671,
                    0.6489481,
                    0.63768298,
                    0.62716971,
                    0.61876397,
                    0.61082755,
                    0.60427772,
                    0.59741976,
                    0.59177741,
                    0.5870773,
                    0.58292462,
                    0.58141689,
                    0.57973476,
                    0.58049471,
                    0.58280279,
                    0.58561492,
                    0.58979099,
                    0.59310853,
                    0.59700109,
                    0.60157219,
                    0.60336097,
                    0.60555331,
                    0.60896068,
                    0.61337866,
                    0.61852465,
                    0.62655929,
                    0.63707128,
                    0.6483534,
                    0.6587092,
                    0.66674618,
                    0.66798851,
                    0.65925168,
                    0.64099533,
                    0.61519263,
                    0.5829609,
                    0.55150764,
                    0.52589593,
                    0.50665129,
                    0.49612167,
                    0.49873702,
                    0.5117356,
                    0.52875232,
                    0.54241942,
                    0.53768022,
                    0.49573105,
                    0.41916397,
                    0.32670548,
                    0.23104246,
                    0.14852103,
                    0.08967661,
                    0.05496955,
                    0.03325212,
                    0.01976446,
                    0.00783771,
                    0.00128398,
                ]
                self.specBandInfo["B8"] = specBandObj

                # B8A
                specBandObj = ARCSISen2SpectralBandObj()
                specBandObj.phyBandName = "B8A"
                specBandObj.bandID = "8"
                specBandObj.imgRes = 20
                specBandObj.wvLenMin = 848
                specBandObj.wvLenMax = 881
                specBandObj.wvLenCen = 864.8
                specBandObj.respFuncStep = 1
                specBandObj.respFunc = [
                    0.0016587,
                    0.01322143,
                    0.02469164,
                    0.05133023,
                    0.10485306,
                    0.21639327,
                    0.38460415,
                    0.58535033,
                    0.77394613,
                    0.87784514,
                    0.91437737,
                    0.92209877,
                    0.92564458,
                    0.9293724,
                    0.93569013,
                    0.94639017,
                    0.95565571,
                    0.96536061,
                    0.97439721,
                    0.97984594,
                    0.98330113,
                    0.98288901,
                    0.98846942,
                    1,
                    0.99957999,
                    0.92089575,
                    0.72838861,
                    0.47188018,
                    0.23786107,
                    0.10682374,
                    0.04603695,
                    0.02219884,
                    0.00879487,
                    0.00046171,
                ]
                self.specBandInfo["B8A"] = specBandObj

                # B09
                specBandObj = ARCSISen2SpectralBandObj()
                specBandObj.phyBandName = "B9"
                specBandObj.bandID = "9"
                specBandObj.imgRes = 60
                specBandObj.wvLenMin = 932
                specBandObj.wvLenMax = 958
                specBandObj.wvLenCen = 945
                specBandObj.respFuncStep = 1
                specBandObj.respFunc = [
                    0.01805614,
                    0.06583501,
                    0.18513673,
                    0.40896107,
                    0.6807859,
                    0.87492845,
                    0.93105831,
                    0.96430107,
                    0.98449689,
                    0.99148444,
                    0.99741262,
                    0.97773458,
                    0.9794157,
                    0.99836495,
                    0.98976032,
                    1,
                    0.98740831,
                    0.98535381,
                    0.95618373,
                    0.96549887,
                    0.93078391,
                    0.86340691,
                    0.70418342,
                    0.44996198,
                    0.20134116,
                    0.05969267,
                    0.0138846,
                ]
                self.specBandInfo["B9"] = specBandObj

                # B10
                specBandObj = ARCSISen2SpectralBandObj()
                specBandObj.phyBandName = "B10"
                specBandObj.bandID = "10"
                specBandObj.imgRes = 60
                specBandObj.wvLenMin = 1337
                specBandObj.wvLenMax = 1412
                specBandObj.wvLenCen = 1373.5
                specBandObj.respFuncStep = 1
                specBandObj.respFunc = [
                    0.00024052,
                    5.404e-05,
                    3.052e-05,
                    2.872e-05,
                    7.632e-05,
                    0.00010949,
                    8.804e-05,
                    0.00012356,
                    0.00017424,
                    0.0003317,
                    0.00036891,
                    0.0004467,
                    0.00065919,
                    0.0010913,
                    0.00196903,
                    0.00373668,
                    0.00801754,
                    0.01884719,
                    0.04466732,
                    0.10165546,
                    0.20111776,
                    0.34284841,
                    0.50710992,
                    0.6632068,
                    0.78377143,
                    0.86153862,
                    0.91000261,
                    0.94193255,
                    0.96182259,
                    0.97365119,
                    0.98169786,
                    0.98795826,
                    0.99283342,
                    0.99649788,
                    0.99906011,
                    1,
                    0.99907734,
                    0.99601604,
                    0.9909083,
                    0.98479854,
                    0.97802142,
                    0.97030114,
                    0.96080954,
                    0.94849765,
                    0.93314108,
                    0.91482336,
                    0.8937997,
                    0.86825426,
                    0.83023193,
                    0.76384193,
                    0.65440009,
                    0.50671604,
                    0.35014737,
                    0.21799972,
                    0.12643091,
                    0.06768988,
                    0.0322709,
                    0.013544,
                    0.00544557,
                    0.00237642,
                    0.00111267,
                    0.00053796,
                    0.0003457,
                    0.00017488,
                    0.00021619,
                    0.00019479,
                    0.00010421,
                    5.919e-05,
                    5.109e-05,
                    6.115e-05,
                    5.527e-05,
                    3.856e-05,
                    3.147e-05,
                    0.00012289,
                    0.0001089,
                    2.502e-05,
                ]
                self.specBandInfo["B10"] = specBandObj

                # B11
                specBandObj = ARCSISen2SpectralBandObj()
                specBandObj.phyBandName = "B11"
                specBandObj.bandID = "11"
                specBandObj.imgRes = 20
                specBandObj.wvLenMin = 1539
                specBandObj.wvLenMax = 1682
                specBandObj.wvLenCen = 1613.7
                specBandObj.respFuncStep = 1
                specBandObj.respFunc = [
                    6.79e-06,
                    6.66e-06,
                    8e-06,
                    2.734e-05,
                    3.685e-05,
                    8.851e-05,
                    0.00014522,
                    0.00024812,
                    0.00047627,
                    0.00056335,
                    0.00065326,
                    0.00089835,
                    0.00114664,
                    0.00165604,
                    0.00241611,
                    0.00350246,
                    0.00524274,
                    0.0081538,
                    0.01237062,
                    0.0186097,
                    0.02721853,
                    0.03879155,
                    0.05379167,
                    0.07353187,
                    0.09932758,
                    0.1334178,
                    0.18029249,
                    0.24484994,
                    0.32834511,
                    0.42749961,
                    0.53576798,
                    0.64570396,
                    0.74245998,
                    0.81447017,
                    0.85866596,
                    0.87924777,
                    0.88665266,
                    0.888727,
                    0.89105732,
                    0.89725046,
                    0.90632982,
                    0.91627527,
                    0.9263751,
                    0.93515828,
                    0.94226446,
                    0.94739906,
                    0.95131987,
                    0.95416808,
                    0.95635128,
                    0.95813297,
                    0.96062738,
                    0.96344083,
                    0.96577764,
                    0.96818134,
                    0.97104025,
                    0.97343195,
                    0.97597444,
                    0.97865413,
                    0.97994672,
                    0.98064126,
                    0.98094979,
                    0.98143338,
                    0.98123856,
                    0.98068083,
                    0.98033995,
                    0.98101894,
                    0.98268503,
                    0.98507875,
                    0.98777658,
                    0.9903608,
                    0.99202087,
                    0.9933069,
                    0.99256744,
                    0.99044883,
                    0.98717314,
                    0.98353656,
                    0.9800432,
                    0.97617287,
                    0.97253451,
                    0.96977033,
                    0.96762556,
                    0.9662626,
                    0.96572411,
                    0.96592079,
                    0.96729798,
                    0.96975438,
                    0.97337748,
                    0.97862858,
                    0.98345358,
                    0.98765317,
                    0.9919238,
                    0.99554959,
                    0.99767411,
                    0.99866451,
                    0.99941783,
                    0.99930984,
                    0.99885298,
                    0.99913515,
                    0.99973164,
                    0.99973592,
                    1,
                    0.9998438,
                    0.9967639,
                    0.99175576,
                    0.9859206,
                    0.97887302,
                    0.97029262,
                    0.96135891,
                    0.95379752,
                    0.94709017,
                    0.94228614,
                    0.93919512,
                    0.93616637,
                    0.92889205,
                    0.9129921,
                    0.88158383,
                    0.82602164,
                    0.74412949,
                    0.64281662,
                    0.53483955,
                    0.42772166,
                    0.32439525,
                    0.23488131,
                    0.16445229,
                    0.11056237,
                    0.07271886,
                    0.04634859,
                    0.02949618,
                    0.01941871,
                    0.0133487,
                    0.00934594,
                    0.00654231,
                    0.00487921,
                    0.00341903,
                    0.00249864,
                    0.00196431,
                    0.00142754,
                    0.00105878,
                    0.00049978,
                    0.00022833,
                    0.00015999,
                    3.415e-05,
                    4.517e-05,
                    1.313e-05,
                ]
                self.specBandInfo["B11"] = specBandObj

                # B12
                specBandObj = ARCSISen2SpectralBandObj()
                specBandObj.phyBandName = "B12"
                specBandObj.bandID = "12"
                specBandObj.imgRes = 20
                specBandObj.wvLenMin = 2078
                specBandObj.wvLenMax = 2320
                specBandObj.wvLenCen = 2202.4
                specBandObj.respFuncStep = 1
                specBandObj.respFunc = [
                    0.00063835,
                    0.00102286,
                    0.00288712,
                    0.00399879,
                    0.00658916,
                    0.00765458,
                    0.00799918,
                    0.00853524,
                    0.00929493,
                    0.00999614,
                    0.01096645,
                    0.01208363,
                    0.01335837,
                    0.01501119,
                    0.01711931,
                    0.01977307,
                    0.02332743,
                    0.02765779,
                    0.03320435,
                    0.04020464,
                    0.04886709,
                    0.0596238,
                    0.07315348,
                    0.09050885,
                    0.11143964,
                    0.13686671,
                    0.16776886,
                    0.20341457,
                    0.24281992,
                    0.28484195,
                    0.32711894,
                    0.36834301,
                    0.40794043,
                    0.4447145,
                    0.47647207,
                    0.50303896,
                    0.52524762,
                    0.54328057,
                    0.55717994,
                    0.5685619,
                    0.57895708,
                    0.58860881,
                    0.59881758,
                    0.60990899,
                    0.62128986,
                    0.63421311,
                    0.64847648,
                    0.66363778,
                    0.67997936,
                    0.69609688,
                    0.71189957,
                    0.7269499,
                    0.74124079,
                    0.75734734,
                    0.77201504,
                    0.78552587,
                    0.79818641,
                    0.80962939,
                    0.81965718,
                    0.82855741,
                    0.83668178,
                    0.84440292,
                    0.85106862,
                    0.85321701,
                    0.85471321,
                    0.8561428,
                    0.85778963,
                    0.8594989,
                    0.86142876,
                    0.86322831,
                    0.86511218,
                    0.8672932,
                    0.86967076,
                    0.87427502,
                    0.87856212,
                    0.88241466,
                    0.88590611,
                    0.8894516,
                    0.89320419,
                    0.8966738,
                    0.89987484,
                    0.90257636,
                    0.90481219,
                    0.90550545,
                    0.90564491,
                    0.90548208,
                    0.90513822,
                    0.90476379,
                    0.90406427,
                    0.90332978,
                    0.90274309,
                    0.90235795,
                    0.90196488,
                    0.90340528,
                    0.90429478,
                    0.90529761,
                    0.90642862,
                    0.90807348,
                    0.91010493,
                    0.91293181,
                    0.91556686,
                    0.91842631,
                    0.92128288,
                    0.92431702,
                    0.92719913,
                    0.92972159,
                    0.93190455,
                    0.93412538,
                    0.93588954,
                    0.93707083,
                    0.93762594,
                    0.93828534,
                    0.93763643,
                    0.94042634,
                    0.94250397,
                    0.94324531,
                    0.94301861,
                    0.94210283,
                    0.94061808,
                    0.93841726,
                    0.93665003,
                    0.93524569,
                    0.93301102,
                    0.92686708,
                    0.92104485,
                    0.91547175,
                    0.91100989,
                    0.90828339,
                    0.9072733,
                    0.90817907,
                    0.91115631,
                    0.91617845,
                    0.92284525,
                    0.92059829,
                    0.91947472,
                    0.91947973,
                    0.92126575,
                    0.92451632,
                    0.92772589,
                    0.93196884,
                    0.93676408,
                    0.94147739,
                    0.94679545,
                    0.95119533,
                    0.95443018,
                    0.95704142,
                    0.95972628,
                    0.9625372,
                    0.96485326,
                    0.96603599,
                    0.96664138,
                    0.96630455,
                    0.96545713,
                    0.96484036,
                    0.96365512,
                    0.96169531,
                    0.95944859,
                    0.95732078,
                    0.95513625,
                    0.95355574,
                    0.95273072,
                    0.95217795,
                    0.95172542,
                    0.9521403,
                    0.95263595,
                    0.95405248,
                    0.95707559,
                    0.96063594,
                    0.96421772,
                    0.96830187,
                    0.97268597,
                    0.97741944,
                    0.98289489,
                    0.9871429,
                    0.99073348,
                    0.99398244,
                    0.99678431,
                    0.99875181,
                    1,
                    0.9999284,
                    0.9991523,
                    0.99712951,
                    0.99388228,
                    0.98968273,
                    0.98373274,
                    0.97621057,
                    0.96780985,
                    0.95833495,
                    0.94842856,
                    0.93818752,
                    0.9277078,
                    0.91702104,
                    0.90597951,
                    0.89384371,
                    0.88165575,
                    0.86861704,
                    0.85460324,
                    0.84058628,
                    0.82598123,
                    0.80948042,
                    0.79182917,
                    0.7724052,
                    0.74907137,
                    0.72031195,
                    0.68815487,
                    0.65125598,
                    0.6100244,
                    0.56600904,
                    0.52095058,
                    0.47464344,
                    0.42924778,
                    0.38584718,
                    0.34208462,
                    0.30067509,
                    0.26317221,
                    0.22770037,
                    0.19571781,
                    0.16808736,
                    0.14467686,
                    0.12482737,
                    0.10823403,
                    0.09439655,
                    0.08235799,
                    0.07149445,
                    0.0626855,
                    0.05498009,
                    0.04818852,
                    0.04285814,
                    0.03859244,
                    0.03494044,
                    0.03199172,
                    0.02958044,
                    0.02741084,
                    0.02556884,
                    0.02395058,
                    0.02166741,
                    0.0191457,
                    0.01632139,
                    0.0109837,
                    0.00736032,
                    0.00649061,
                    0.00469736,
                    0.00205874,
                ]
                self.specBandInfo["B12"] = specBandObj

            tree = None

            ####### READ GRANULE HEADER FILES ##########
            granuleHdr = ""
            if (hdrFileVersion == "psd14") and (not mixedPSD14PSD13):
                granuleDIR = os.path.join(self.sen2FileBaseDIR, "GRANULE")
                granLC1DIRsIn = os.listdir(granuleDIR)
                granLC1DIRs = []
                for dirStr in granLC1DIRsIn:
                    if "L1C" in dirStr:
                        granLC1DIRs.append(dirStr)
                if len(granLC1DIRs) != 1:
                    raise ARCSIException("Couldn't find the granule directory")
                granuleHdr = os.path.join(
                    os.path.join(granuleDIR, granLC1DIRs[0]), "MTD_TL.xml"
                )
            elif (hdrFileVersion == "psd13") or (
                mixedPSD14PSD13 and (hdrFileVersion == "psd14")
            ):
                granuleDIR = os.path.join(
                    os.path.join(self.sen2FileBaseDIR, "GRANULE"), granuleIdentifier
                )
                granuleTmpFiles = glob.glob(
                    os.path.join(granuleDIR, "S2A_OPER_MTD_L1C*.xml")
                )
                if len(granuleTmpFiles) == 1:
                    granuleHdr = granuleTmpFiles[0]
                else:
                    raise ARCSIException("Cannot find the granule header file.")
            else:
                raise ARCSIException(
                    "Do not recognised header format: '" + hdrFileVersion + "'"
                )

            tree = ET.parse(granuleHdr)
            root = tree.getroot()

            # Get Geometric Info Tag
            geometricInfoTag = None

            # Try for diffrent versions
            for hdr_format_version in [12, 14, 15]:
                if geometricInfoTag is None:
                    geometricInfoTag = root.find(
                        f"{{https://psd-{hdr_format_version}.sentinel2.eo.esa.int/PSD/S2_PDI_Level-1C_Tile_Metadata.xsd}}Geometric_Info"
                    )

            # If not found raise exception
            if geometricInfoTag is None:
                raise ARCSIException(
                    "Geometric_Info in granular header filed does not match expected formats"
                )

            # Get Tile Geocoding tag.
            tileGeocoding = geometricInfoTag.find("Tile_Geocoding")

            # Get EPSG Code.
            epsgCodeStr = tileGeocoding.find("HORIZONTAL_CS_CODE").text.strip()
            self.epsgCode = rsgislib.tools.utils.str_to_int(epsgCodeStr.split(":")[1])
            inProj = osr.SpatialReference()
            inProj.ImportFromEPSG(self.epsgCode)
            if self.inWKT == "":
                self.inWKT = inProj.ExportToWkt()

            self.projNameStr = ""
            utmZone = inProj.GetUTMZone()
            utmHem = ""
            if utmZone < 0:
                utmHem = "s"
            else:
                utmHem = "n"
            if utmZone == 0:
                raise ARCSIException(
                    "Sen2 Image is not projected with UTM - contact support as header not currently supported."
                )
            else:
                utmZoneStr = str(utmZone).replace("-", "")
                self.projNameStr = "utm" + utmZoneStr + utmHem

            tlXRes10 = 0.0
            tlYRes10 = 0.0
            xSizeRes10 = 0
            ySizeRes10 = 0
            tlXRes20 = 0.0
            tlYRes20 = 0.0
            xSizeRes20 = 0
            ySizeRes20 = 0
            tlXRes60 = 0.0
            tlYRes60 = 0.0
            xSizeRes60 = 0
            ySizeRes60 = 0
            for tileGeocodeChild in tileGeocoding:
                if tileGeocodeChild.tag == "Geoposition":
                    if tileGeocodeChild.attrib["resolution"] == "10":
                        tlXRes10 = rsgislib.tools.utils.str_to_float(
                            tileGeocodeChild.find("ULX").text.strip()
                        )
                        tlYRes10 = rsgislib.tools.utils.str_to_float(
                            tileGeocodeChild.find("ULY").text.strip()
                        )
                    elif tileGeocodeChild.attrib["resolution"] == "20":
                        tlXRes20 = rsgislib.tools.utils.str_to_float(
                            tileGeocodeChild.find("ULX").text.strip()
                        )
                        tlYRes20 = rsgislib.tools.utils.str_to_float(
                            tileGeocodeChild.find("ULY").text.strip()
                        )
                    elif tileGeocodeChild.attrib["resolution"] == "60":
                        tlXRes60 = rsgislib.tools.utils.str_to_float(
                            tileGeocodeChild.find("ULX").text.strip()
                        )
                        tlYRes60 = rsgislib.tools.utils.str_to_float(
                            tileGeocodeChild.find("ULY").text.strip()
                        )
                    else:
                        raise ARCSIException("Resolution must be 10, 20 or 60m")
                elif tileGeocodeChild.tag == "Size":
                    if tileGeocodeChild.attrib["resolution"] == "10":
                        xSizeRes10 = rsgislib.tools.utils.str_to_int(
                            tileGeocodeChild.find("NCOLS").text.strip()
                        )
                        ySizeRes10 = rsgislib.tools.utils.str_to_int(
                            tileGeocodeChild.find("NROWS").text.strip()
                        )
                    elif tileGeocodeChild.attrib["resolution"] == "20":
                        xSizeRes20 = rsgislib.tools.utils.str_to_int(
                            tileGeocodeChild.find("NCOLS").text.strip()
                        )
                        ySizeRes20 = rsgislib.tools.utils.str_to_int(
                            tileGeocodeChild.find("NROWS").text.strip()
                        )
                    elif tileGeocodeChild.attrib["resolution"] == "60":
                        xSizeRes60 = rsgislib.tools.utils.str_to_int(
                            tileGeocodeChild.find("NCOLS").text.strip()
                        )
                        ySizeRes60 = rsgislib.tools.utils.str_to_int(
                            tileGeocodeChild.find("NROWS").text.strip()
                        )
                    else:
                        raise ARCSIException("Resolution must be 10, 20 or 60m")

            brXRes10 = tlXRes10 + (xSizeRes10 * 10)
            brYRes10 = tlYRes10 - (ySizeRes10 * 10)
            brXRes20 = tlXRes20 + (xSizeRes20 * 20)
            brYRes20 = tlYRes20 - (ySizeRes20 * 20)
            brXRes60 = tlXRes60 + (xSizeRes60 * 60)
            brYRes60 = tlYRes60 - (ySizeRes60 * 60)

            self.xTL = min([tlXRes10, tlXRes20, tlXRes60])
            self.yTL = max([tlYRes10, tlYRes20, tlYRes60])
            self.xBR = max([brXRes10, brXRes20, brXRes60])
            self.yBR = min([brYRes10, brYRes20, brYRes60])
            self.xTR = self.xBR
            self.yTR = self.yTL
            self.xBL = self.xTL
            self.yBL = self.yBR

            self.xCentre = self.xTL + ((self.xTR - self.xTL) / 2)
            self.yCentre = self.yBR + ((self.yTL - self.yBR) / 2)

            (
                self.lonCentre,
                self.latCentre,
            ) = rsgislib.tools.geometrytools.reproj_point_to_wgs84(
                inProj, self.xCentre, self.yCentre
            )
            self.lonTL, self.latTL = rsgislib.tools.geometrytools.reproj_point_to_wgs84(
                inProj, self.xTL, self.yTL
            )
            self.lonTR, self.latTR = rsgislib.tools.geometrytools.reproj_point_to_wgs84(
                inProj, self.xTR, self.yTR
            )
            self.lonBL, self.latBL = rsgislib.tools.geometrytools.reproj_point_to_wgs84(
                inProj, self.xBL, self.yBL
            )
            self.lonBR, self.latBR = rsgislib.tools.geometrytools.reproj_point_to_wgs84(
                inProj, self.xBR, self.yBR
            )

            # Get Tile angles tag.
            tileAngles = geometricInfoTag.find("Tile_Angles")

            sunAngleTag = tileAngles.find("Mean_Sun_Angle")
            self.solarZenith = rsgislib.tools.utils.str_to_float(
                sunAngleTag.find("ZENITH_ANGLE").text.strip()
            )
            self.solarAzimuth = rsgislib.tools.utils.str_to_float(
                sunAngleTag.find("AZIMUTH_ANGLE").text.strip()
            )

            senZenVals = []
            senAzVals = []
            for meanViewIncAngleTag in tileAngles.find(
                "Mean_Viewing_Incidence_Angle_List"
            ):
                if meanViewIncAngleTag.tag == "Mean_Viewing_Incidence_Angle":
                    senZenVals.append(
                        rsgislib.tools.utils.str_to_float(
                            meanViewIncAngleTag.find("ZENITH_ANGLE").text.strip()
                        )
                    )
                    senAzVals.append(
                        rsgislib.tools.utils.str_to_float(
                            meanViewIncAngleTag.find("AZIMUTH_ANGLE").text.strip()
                        )
                    )

            self.sensorZenith = statistics.mean(senZenVals)
            self.sensorAzimuth = statistics.mean(senAzVals)

            ################################
            #### CODE FROM PYTHON FMASK ####
            ################################
            # Upper-left corners of images at different resolutions. As far as I can
            # work out, these coords appear to be the upper left corner of the upper left
            # pixel, i.e. equivalent to GDAL's convention. This also means that they
            # are the same for the different resolutions, which is nice.
            self.ulxyByRes = {}
            posNodeList = tileGeocoding.findall("Geoposition")
            for posNode in posNodeList:
                res = posNode.attrib["resolution"]
                ulx = float(posNode.find("ULX").text)
                uly = float(posNode.find("ULY").text)
                self.ulxyByRes[res] = (ulx, uly)

            # Sun and satellite angles.
            sunZenithNode = tileAngles.find("Sun_Angles_Grid").find("Zenith")
            self.angleGridXres = float(sunZenithNode.find("COL_STEP").text)
            self.angleGridYres = float(sunZenithNode.find("ROW_STEP").text)
            self.sunZenithGrid = self.makeValueArray(sunZenithNode.find("Values_List"))
            sunAzimuthNode = tileAngles.find("Sun_Angles_Grid").find("Azimuth")
            self.sunAzimuthGrid = self.makeValueArray(
                sunAzimuthNode.find("Values_List")
            )
            self.anglesGridShape = self.sunAzimuthGrid.shape

            # Now build up the viewing angle per grid cell, from the separate layers
            # given for each detector for each band. Initially I am going to keep
            # the bands separate, just to see how that looks.
            # The names of things in the XML suggest that these are view angles,
            # but the numbers suggest that they are angles as seen from the pixel's
            # frame of reference on the ground, i.e. they are in fact what we ultimately want.
            viewingAngleNodeList = tileAngles.findall("Viewing_Incidence_Angles_Grids")
            self.viewZenithDict = self.buildViewAngleArr(viewingAngleNodeList, "Zenith")
            self.viewAzimuthDict = self.buildViewAngleArr(
                viewingAngleNodeList, "Azimuth"
            )

            # Make a guess at the coordinates of the angle grids. These are not given
            # explicitly in the XML, and don't line up exactly with the other grids, so I am
            # making a rough estimate. Because the angles don't change rapidly across these
            # distances, it is not important if I am a bit wrong (although it would be nice
            # to be exactly correct!).
            (ulx, uly) = self.ulxyByRes["10"]
            self.anglesULXY = (
                ulx - self.angleGridXres / 2.0,
                uly + self.angleGridYres / 2.0,
            )
            ################################

            # Get 6S spectral response functions.
            for specRspBand in self.specBandInfo:
                respFuncSize = len(self.specBandInfo[specRspBand].respFunc)
                tmpWvLenMin = self.specBandInfo[specRspBand].wvLenMin
                tmpWvLenMax = self.specBandInfo[specRspBand].wvLenMax
                numWvLens = math.ceil(tmpWvLenMax - tmpWvLenMin)
                lenDiff = respFuncSize - numWvLens
                if lenDiff > 0:
                    tmpWvLenMax = tmpWvLenMax + lenDiff
                elif lenDiff < 0:
                    tmpWvLenMin = tmpWvLenMin + (lenDiff * (-1))
                wvLensIn = numpy.arange(tmpWvLenMin, tmpWvLenMax, 1)
                olWVLens, olSpecResp = arcsilib.arcsiutils.resampleSpectralResponseFunc(
                    wvLensIn, self.specBandInfo[specRspBand].respFunc, 2.5, "linear"
                )
                self.specBandInfo[specRspBand].respFunc6S = olSpecResp
                self.specBandInfo[specRspBand].wvLenMin6S = olWVLens[0] / 1000
                self.specBandInfo[specRspBand].wvLenMax6S = olWVLens[-1] / 1000

        except Exception as e:
            raise e

    def getSolarIrrStdSolarGeom(self):
        """
        Get Solar Azimuth and Zenith as standard geometry.
        Azimuth: N=0, E=90, S=180, W=270.
        """
        return (self.solarAzimuth, self.solarZenith)

    def getSensorViewGeom(self):
        """
        Get sensor viewing angles
        returns (viewAzimuth, viewZenith)
        """
        return (0.0, 0.0)

    def generateOutputBaseName(self):
        """
        Customises the generic name for the WorldView2 sensor
        """
        outname = self.defaultGenBaseOutFileName()
        orbNumStr = str(self.orbitNumber)
        if len(orbNumStr) == 1:
            orbNumStr = "00" + orbNumStr
        if len(orbNumStr) == 2:
            orbNumStr = "0" + orbNumStr
        gen_time_str = self.generationTime.strftime("%Y%m%d%H%M%S")
        outname = (
            outname
            + "_"
            + self.uniqueTileID
            + "_ORB"
            + orbNumStr
            + "_"
            + gen_time_str
            + "_"
            + self.projNameStr
        )
        return outname

    def expectedImageDataPresent(self):
        imageDataPresent = True

        if not os.path.exists(self.sen2ImgB01):
            imageDataPresent = False
        elif not os.path.exists(self.sen2ImgB02):
            imageDataPresent = False
        elif not os.path.exists(self.sen2ImgB03):
            imageDataPresent = False
        elif not os.path.exists(self.sen2ImgB04):
            imageDataPresent = False
        elif not os.path.exists(self.sen2ImgB05):
            imageDataPresent = False
        elif not os.path.exists(self.sen2ImgB06):
            imageDataPresent = False
        elif not os.path.exists(self.sen2ImgB07):
            imageDataPresent = False
        elif not os.path.exists(self.sen2ImgB8A):
            imageDataPresent = False
        elif not os.path.exists(self.sen2ImgB08):
            imageDataPresent = False
        elif not os.path.exists(self.sen2ImgB09):
            imageDataPresent = False
        elif not os.path.exists(self.sen2ImgB10):
            imageDataPresent = False
        elif not os.path.exists(self.sen2ImgB11):
            imageDataPresent = False
        elif not os.path.exists(self.sen2ImgB12):
            imageDataPresent = False

        return imageDataPresent

    def imgNeedMosaicking(self):
        return False

    def inImgsDiffRes(self):
        return True

    def applyImageDataMask(
        self,
        inputHeader,
        inputImage,
        outputPath,
        outputMaskName,
        outputImgName,
        outFormat,
        outWKTFile,
    ):
        raise ARCSIException(
            "Sentinel-2 does not provide any image masks, do not use the MASK option."
        )

    def mosaicImageTiles(self, outputPath):
        raise ARCSIException("Image data does not need mosaicking")

    def resampleImgRes(
        self,
        outputPath,
        resampleToLowResImg,
        resampleMethod=rsgislib.INTERP_CUBIC,
        multicore=False,
    ):
        outBaseName = self.generateOutputBaseName()
        if resampleToLowResImg:
            # Resample to 20 m
            self.sen2ImgB02_20m = os.path.join(outputPath, outBaseName + "_B02_20m.kea")
            rsgislib.imageutils.resample_img_to_match(
                self.sen2ImgB05,
                self.sen2ImgB02,
                self.sen2ImgB02_20m,
                "KEA",
                resampleMethod,
                rsgislib.TYPE_16UINT,
                0.0,
                multicore,
            )
            self.sen2ImgB03_20m = os.path.join(outputPath, outBaseName + "_B03_20m.kea")
            rsgislib.imageutils.resample_img_to_match(
                self.sen2ImgB05,
                self.sen2ImgB03,
                self.sen2ImgB03_20m,
                "KEA",
                resampleMethod,
                rsgislib.TYPE_16UINT,
                0.0,
                multicore,
            )
            self.sen2ImgB04_20m = os.path.join(outputPath, outBaseName + "_B04_20m.kea")
            rsgislib.imageutils.resample_img_to_match(
                self.sen2ImgB05,
                self.sen2ImgB04,
                self.sen2ImgB04_20m,
                "KEA",
                resampleMethod,
                rsgislib.TYPE_16UINT,
                0.0,
                multicore,
            )
            self.sen2ImgB08_20m = os.path.join(outputPath, outBaseName + "_B08_20m.kea")
            rsgislib.imageutils.resample_img_to_match(
                self.sen2ImgB05,
                self.sen2ImgB08,
                self.sen2ImgB08_20m,
                "KEA",
                resampleMethod,
                rsgislib.TYPE_16UINT,
                0.0,
                multicore,
            )
            self.resampleTo20m = True
        else:
            # Resample to 10 m
            self.sen2ImgB05_10m = os.path.join(outputPath, outBaseName + "_B05_10m.kea")
            rsgislib.imageutils.resample_img_to_match(
                self.sen2ImgB02,
                self.sen2ImgB05,
                self.sen2ImgB05_10m,
                "KEA",
                resampleMethod,
                rsgislib.TYPE_16UINT,
                0.0,
                multicore,
            )
            self.sen2ImgB06_10m = os.path.join(outputPath, outBaseName + "_B06_10m.kea")
            rsgislib.imageutils.resample_img_to_match(
                self.sen2ImgB02,
                self.sen2ImgB06,
                self.sen2ImgB06_10m,
                "KEA",
                resampleMethod,
                rsgislib.TYPE_16UINT,
                0.0,
                multicore,
            )
            self.sen2ImgB07_10m = os.path.join(outputPath, outBaseName + "_B07_10m.kea")
            rsgislib.imageutils.resample_img_to_match(
                self.sen2ImgB02,
                self.sen2ImgB07,
                self.sen2ImgB07_10m,
                "KEA",
                resampleMethod,
                rsgislib.TYPE_16UINT,
                0.0,
                multicore,
            )
            self.sen2ImgB8A_10m = os.path.join(
                outputPath, outBaseName + "_B08A_10m.kea"
            )
            rsgislib.imageutils.resample_img_to_match(
                self.sen2ImgB02,
                self.sen2ImgB8A,
                self.sen2ImgB8A_10m,
                "KEA",
                resampleMethod,
                rsgislib.TYPE_16UINT,
                0.0,
                multicore,
            )
            self.sen2ImgB11_10m = os.path.join(outputPath, outBaseName + "_B11_10m.kea")
            rsgislib.imageutils.resample_img_to_match(
                self.sen2ImgB02,
                self.sen2ImgB11,
                self.sen2ImgB11_10m,
                "KEA",
                resampleMethod,
                rsgislib.TYPE_16UINT,
                0.0,
                multicore,
            )
            self.sen2ImgB12_10m = os.path.join(outputPath, outBaseName + "_B12_10m.kea")
            rsgislib.imageutils.resample_img_to_match(
                self.sen2ImgB02,
                self.sen2ImgB12,
                self.sen2ImgB12_10m,
                "KEA",
                resampleMethod,
                rsgislib.TYPE_16UINT,
                0.0,
                multicore,
            )
            self.resampleTo20m = False

    def sharpenLowResRadImgBands(self, inputImg, outputImage, outFormat):
        bandInfo = []
        bandInfo.append(
            rsgislib.imageutils.SharpBandInfo(
                band=1, status=rsgislib.SHARP_RES_HIGH, name="Blue"
            )
        )
        bandInfo.append(
            rsgislib.imageutils.SharpBandInfo(
                band=2, status=rsgislib.SHARP_RES_HIGH, name="Green"
            )
        )
        bandInfo.append(
            rsgislib.imageutils.SharpBandInfo(
                band=3, status=rsgislib.SHARP_RES_HIGH, name="Red"
            )
        )
        bandInfo.append(
            rsgislib.imageutils.SharpBandInfo(
                band=4, status=rsgislib.SHARP_RES_LOW, name="RE_B5"
            )
        )
        bandInfo.append(
            rsgislib.imageutils.SharpBandInfo(
                band=5, status=rsgislib.SHARP_RES_LOW, name="RE_B6"
            )
        )
        bandInfo.append(
            rsgislib.imageutils.SharpBandInfo(
                band=6, status=rsgislib.SHARP_RES_LOW, name="RE_B7"
            )
        )
        bandInfo.append(
            rsgislib.imageutils.SharpBandInfo(
                band=7, status=rsgislib.SHARP_RES_HIGH, name="NIR_B8"
            )
        )
        bandInfo.append(
            rsgislib.imageutils.SharpBandInfo(
                band=8, status=rsgislib.SHARP_RES_LOW, name="NIR_B8A"
            )
        )
        bandInfo.append(
            rsgislib.imageutils.SharpBandInfo(
                band=9, status=rsgislib.SHARP_RES_LOW, name="SWIR1"
            )
        )
        bandInfo.append(
            rsgislib.imageutils.SharpBandInfo(
                band=10, status=rsgislib.SHARP_RES_LOW, name="SWIR2"
            )
        )
        rsgislib.imageutils.sharpen_low_res_bands(
            inputImg, outputImage, bandInfo, 7, 0, outFormat, rsgislib.TYPE_32FLOAT
        )

    def generateValidImageDataMask(
        self, outputPath, outputMaskName, viewAngleImg, outFormat
    ):
        print("Generate valid image mask")
        # Generate the valid image mask
        outputImage = os.path.join(outputPath, outputMaskName)
        inImgBands = []
        if self.resampleTo20m:
            inImgBands.append(self.sen2ImgB02_20m)
            inImgBands.append(self.sen2ImgB03_20m)
            inImgBands.append(self.sen2ImgB04_20m)
            inImgBands.append(self.sen2ImgB05)
            inImgBands.append(self.sen2ImgB06)
            inImgBands.append(self.sen2ImgB07)
            inImgBands.append(self.sen2ImgB08_20m)
            inImgBands.append(self.sen2ImgB8A)
            inImgBands.append(self.sen2ImgB11)
            inImgBands.append(self.sen2ImgB12)
        else:
            inImgBands.append(self.sen2ImgB02)
            inImgBands.append(self.sen2ImgB03)
            inImgBands.append(self.sen2ImgB04)
            inImgBands.append(self.sen2ImgB05_10m)
            inImgBands.append(self.sen2ImgB06_10m)
            inImgBands.append(self.sen2ImgB07_10m)
            inImgBands.append(self.sen2ImgB08)
            inImgBands.append(self.sen2ImgB8A_10m)
            inImgBands.append(self.sen2ImgB11_10m)
            inImgBands.append(self.sen2ImgB12_10m)
        rsgislib.imageutils.gen_valid_mask(
            input_imgs=inImgBands,
            output_img=outputImage,
            gdalformat=outFormat,
            no_data_val=self.inNoDataVal,
        )

        ################################
        #### CODE FROM PYTHON FMASK ####
        ################################
        tmpViewAngleImg = os.path.join(
            outputPath, os.path.splitext(os.path.basename(viewAngleImg))[0] + "tmp.kea"
        )
        # Generate the acquasition anagles image.
        # This scale value will convert between DN and radians in output
        # image file, radians = dn * SCALE_TO_RADIANS
        SCALE_TO_RADIANS = 0.01

        drvr = gdal.GetDriverByName("KEA")
        (nrows, ncols) = self.anglesGridShape
        ds = drvr.Create(tmpViewAngleImg, ncols, nrows, 4, gdal.GDT_Int16)
        gt = (
            self.anglesULXY[0],
            self.angleGridXres,
            0,
            self.anglesULXY[1],
            0.0,
            -self.angleGridYres,
        )
        ds.SetGeoTransform(gt)
        sr = osr.SpatialReference()
        sr.ImportFromEPSG(int(self.epsgCode))
        ds.SetProjection(sr.ExportToWkt())

        nullValDN = 1000

        # Get a sorted list of the Sentinel-2 band names. Note that sometimes this
        # is an incomplete list of band names, which appears to be due to a bug in
        # earlier versions of ESA's processing software. I suspect it relates to
        # Anomaly number 11 in the following page.
        # https://sentinel.esa.int/web/sentinel/news/-/article/new-processing-baseline-for-sentinel-2-products
        bandNames = sorted(self.viewAzimuthDict.keys())

        # Mean over all bands
        satAzDeg = numpy.array([self.viewAzimuthDict[i] for i in bandNames])
        satAzDegMeanOverBands = satAzDeg.mean(axis=0)

        satZenDeg = numpy.array([self.viewZenithDict[i] for i in bandNames])
        satZenDegMeanOverBands = satZenDeg.mean(axis=0)

        sunAzDeg = self.sunAzimuthGrid
        sunZenDeg = self.sunZenithGrid

        stackDeg = numpy.array(
            [satAzDegMeanOverBands, satZenDegMeanOverBands, sunAzDeg, sunZenDeg]
        )
        stackRadians = numpy.radians(stackDeg)

        nullmask = numpy.isnan(stackDeg)
        stackRadians[nullmask] = 0.0
        stackDN = numpy.round(stackRadians / SCALE_TO_RADIANS).astype(numpy.int16)
        stackDN[nullmask] = nullValDN

        lnames = ["SatelliteAzimuth", "SatelliteZenith", "SunAzimuth", "SunZenith"]
        for i in range(ds.RasterCount):
            b = ds.GetRasterBand(i + 1)
            b.WriteArray(stackDN[i])
            b.SetNoDataValue(nullValDN)
            b.SetDescription(lnames[i])
        del ds
        ################################
        rsgislib.imageutils.pop_img_stats(tmpViewAngleImg, True, 1000, False)
        rsgislib.imageutils.resample_img_to_match(
            outputImage,
            tmpViewAngleImg,
            viewAngleImg,
            outFormat,
            rsgislib.INTERP_NEAREST_NEIGHBOUR,
            datatype=None,
            multicore=False,
        )
        dataset = gdal.Open(viewAngleImg, gdal.GA_Update)
        if not dataset is None:
            dataset.GetRasterBand(1).SetDescription("SatelliteAzimuth")
            dataset.GetRasterBand(2).SetDescription("SatelliteZenith")
            dataset.GetRasterBand(3).SetDescription("SolorAzimuth")
            dataset.GetRasterBand(4).SetDescription("SolorZenith")
        dataset = None
        drvr.Delete(tmpViewAngleImg)

        return outputImage

    def generateImageSaturationMask(self, outputPath, outputName, outFormat):
        print("Generate Saturation Image")
        outputImage = os.path.join(outputPath, outputName)

        s2Band = collections.namedtuple(
            "S2Band", ["band_name", "input_img", "img_band", "sat_val"]
        )
        bandDefnSeq = list()
        self.inSatDataVal = float(self.inSatDataVal)
        if self.resampleTo20m:
            bandDefnSeq.append(
                s2Band(
                    band_name="Blue",
                    input_img=self.sen2ImgB02_20m,
                    img_band=1,
                    sat_val=self.inSatDataVal,
                )
            )
            bandDefnSeq.append(
                s2Band(
                    band_name="Green",
                    input_img=self.sen2ImgB03_20m,
                    img_band=1,
                    sat_val=self.inSatDataVal,
                )
            )
            bandDefnSeq.append(
                s2Band(
                    band_name="Red",
                    input_img=self.sen2ImgB04_20m,
                    img_band=1,
                    sat_val=self.inSatDataVal,
                )
            )
            bandDefnSeq.append(
                s2Band(
                    band_name="RE_B05",
                    input_img=self.sen2ImgB05,
                    img_band=1,
                    sat_val=self.inSatDataVal,
                )
            )
            bandDefnSeq.append(
                s2Band(
                    band_name="RE_B06",
                    input_img=self.sen2ImgB06,
                    img_band=1,
                    sat_val=self.inSatDataVal,
                )
            )
            bandDefnSeq.append(
                s2Band(
                    band_name="RE_B07",
                    input_img=self.sen2ImgB07,
                    img_band=1,
                    sat_val=self.inSatDataVal,
                )
            )
            bandDefnSeq.append(
                s2Band(
                    band_name="NIR_B08",
                    input_img=self.sen2ImgB08_20m,
                    img_band=1,
                    sat_val=self.inSatDataVal,
                )
            )
            bandDefnSeq.append(
                s2Band(
                    band_name="NIR_B08A",
                    input_img=self.sen2ImgB8A,
                    img_band=1,
                    sat_val=self.inSatDataVal,
                )
            )
            bandDefnSeq.append(
                s2Band(
                    band_name="SWIR1",
                    input_img=self.sen2ImgB11,
                    img_band=1,
                    sat_val=self.inSatDataVal,
                )
            )
            bandDefnSeq.append(
                s2Band(
                    band_name="SWIR2",
                    input_img=self.sen2ImgB12,
                    img_band=1,
                    sat_val=self.inSatDataVal,
                )
            )
        else:
            bandDefnSeq.append(
                s2Band(
                    band_name="Blue",
                    input_img=self.sen2ImgB02,
                    img_band=1,
                    sat_val=self.inSatDataVal,
                )
            )
            bandDefnSeq.append(
                s2Band(
                    band_name="Green",
                    input_img=self.sen2ImgB03,
                    img_band=1,
                    sat_val=self.inSatDataVal,
                )
            )
            bandDefnSeq.append(
                s2Band(
                    band_name="Red",
                    input_img=self.sen2ImgB04,
                    img_band=1,
                    sat_val=self.inSatDataVal,
                )
            )
            bandDefnSeq.append(
                s2Band(
                    band_name="RE_B05",
                    input_img=self.sen2ImgB05_10m,
                    img_band=1,
                    sat_val=self.inSatDataVal,
                )
            )
            bandDefnSeq.append(
                s2Band(
                    band_name="RE_B06",
                    input_img=self.sen2ImgB06_10m,
                    img_band=1,
                    sat_val=self.inSatDataVal,
                )
            )
            bandDefnSeq.append(
                s2Band(
                    band_name="RE_B07",
                    input_img=self.sen2ImgB07_10m,
                    img_band=1,
                    sat_val=self.inSatDataVal,
                )
            )
            bandDefnSeq.append(
                s2Band(
                    band_name="NIR_B08",
                    input_img=self.sen2ImgB08,
                    img_band=1,
                    sat_val=self.inSatDataVal,
                )
            )
            bandDefnSeq.append(
                s2Band(
                    band_name="NIR_B08A",
                    input_img=self.sen2ImgB8A_10m,
                    img_band=1,
                    sat_val=self.inSatDataVal,
                )
            )
            bandDefnSeq.append(
                s2Band(
                    band_name="SWIR1",
                    input_img=self.sen2ImgB11_10m,
                    img_band=1,
                    sat_val=self.inSatDataVal,
                )
            )
            bandDefnSeq.append(
                s2Band(
                    band_name="SWIR2",
                    input_img=self.sen2ImgB12_10m,
                    img_band=1,
                    sat_val=self.inSatDataVal,
                )
            )

        rsgislib.imagecalibration.saturated_pixels_mask(
            outputImage, outFormat, bandDefnSeq
        )
        return outputImage

    def convertImageToRadiance(
        self, outputPath, outputReflName, outputThermalName, outFormat
    ):
        print("Converting to Radiance")
        outputReflImage = os.path.join(outputPath, outputReflName)
        outputThermalImage = None

        img_band02 = self.sen2ImgB02
        img_band03 = self.sen2ImgB03
        img_band04 = self.sen2ImgB04
        img_band05 = self.sen2ImgB05_10m
        img_band06 = self.sen2ImgB06_10m
        img_band07 = self.sen2ImgB07_10m
        img_band08 = self.sen2ImgB08
        img_band08A = self.sen2ImgB8A_10m
        img_band11 = self.sen2ImgB11_10m
        img_band12 = self.sen2ImgB12_10m

        if self.resampleTo20m:
            img_band02 = self.sen2ImgB02_20m
            img_band03 = self.sen2ImgB03_20m
            img_band04 = self.sen2ImgB04_20m
            img_band05 = self.sen2ImgB05
            img_band06 = self.sen2ImgB06
            img_band07 = self.sen2ImgB07
            img_band08 = self.sen2ImgB08_20m
            img_band08A = self.sen2ImgB8A
            img_band11 = self.sen2ImgB11
            img_band12 = self.sen2ImgB12

        img_band02_tmp = None
        img_band03_tmp = None
        img_band04_tmp = None
        img_band05_tmp = None
        img_band06_tmp = None
        img_band07_tmp = None
        img_band08_tmp = None
        img_band08A_tmp = None
        img_band11_tmp = None
        img_band12_tmp = None

        basename = self.generateOutputBaseName() # use a basename which has the generation time for uniqueness

        if abs(self.ratiometric_offs_B1) > 0:
            img_band02_tmp = os.path.join(
                outputPath, f"{basename}_B02_tmpoffapply.kea"
            )
            exp = f"b1 + {self.ratiometric_offs_B1}"
            rsgislib.imagecalc.image_math(
                img_band02, img_band02_tmp, exp, "KEA", rsgislib.TYPE_32FLOAT
            )
            img_band02 = img_band02_tmp

        if abs(self.ratiometric_offs_B2) > 0:
            img_band03_tmp = os.path.join(
                outputPath, f"{basename}_B03_tmpoffapply.kea"
            )
            exp = f"b1 + {self.ratiometric_offs_B2}"
            rsgislib.imagecalc.image_math(
                img_band03, img_band03_tmp, exp, "KEA", rsgislib.TYPE_32FLOAT
            )
            img_band03 = img_band03_tmp

        if abs(self.ratiometric_offs_B3) > 0:
            img_band04_tmp = os.path.join(
                outputPath, f"{basename}_B04_tmpoffapply.kea"
            )
            exp = f"b1 + {self.ratiometric_offs_B3}"
            rsgislib.imagecalc.image_math(
                img_band04, img_band04_tmp, exp, "KEA", rsgislib.TYPE_32FLOAT
            )
            img_band04 = img_band04_tmp

        if abs(self.ratiometric_offs_B4) > 0:
            img_band05_tmp = os.path.join(
                outputPath, f"{basename}_B05_tmpoffapply.kea"
            )
            exp = f"b1 + {self.ratiometric_offs_B4}"
            rsgislib.imagecalc.image_math(
                img_band05, img_band05_tmp, exp, "KEA", rsgislib.TYPE_32FLOAT
            )
            img_band05 = img_band05_tmp

        if abs(self.ratiometric_offs_B5) > 0:
            img_band06_tmp = os.path.join(
                outputPath, f"{basename}_B06_tmpoffapply.kea"
            )
            exp = f"b1 + {self.ratiometric_offs_B5}"
            rsgislib.imagecalc.image_math(
                img_band06, img_band06_tmp, exp, "KEA", rsgislib.TYPE_32FLOAT
            )
            img_band06 = img_band06_tmp

        if abs(self.ratiometric_offs_B6) > 0:
            img_band07_tmp = os.path.join(
                outputPath, f"{basename}_B07_tmpoffapply.kea"
            )
            exp = f"b1 + {self.ratiometric_offs_B6}"
            rsgislib.imagecalc.image_math(
                img_band07, img_band07_tmp, exp, "KEA", rsgislib.TYPE_32FLOAT
            )
            img_band07 = img_band07_tmp

        if abs(self.ratiometric_offs_B7) > 0:
            img_band08_tmp = os.path.join(
                outputPath, f"{basename}_B08_tmpoffapply.kea"
            )
            exp = f"b1 + {self.ratiometric_offs_B7}"
            rsgislib.imagecalc.image_math(
                img_band08, img_band08_tmp, exp, "KEA", rsgislib.TYPE_32FLOAT
            )
            img_band08 = img_band08_tmp

        if abs(self.ratiometric_offs_B8) > 0:
            img_band08A_tmp = os.path.join(
                outputPath, f"{basename}_B08A_tmpoffapply.kea"
            )
            exp = f"b1 + {self.ratiometric_offs_B8}"
            rsgislib.imagecalc.image_math(
                img_band08A, img_band08A_tmp, exp, "KEA", rsgislib.TYPE_32FLOAT
            )
            img_band08A = img_band08A_tmp

        if abs(self.ratiometric_offs_B11) > 0:
            img_band11_tmp = os.path.join(
                outputPath, f"{basename}_B11_tmpoffapply.kea"
            )
            exp = f"b1 + {self.ratiometric_offs_B11}"
            rsgislib.imagecalc.image_math(
                img_band11, img_band11_tmp, exp, "KEA", rsgislib.TYPE_32FLOAT
            )
            img_band11 = img_band11_tmp

        if abs(self.ratiometric_offs_B12) > 0:
            img_band12_tmp = os.path.join(
                outputPath, f"{basename}_B12_tmpoffapply.kea"
            )
            exp = f"b1 + {self.ratiometric_offs_B12}"
            rsgislib.imagecalc.image_math(
                img_band12, img_band12_tmp, exp, "KEA", rsgislib.TYPE_32FLOAT
            )
            img_band12 = img_band12_tmp

        inImgBands = list()
        inImgBands.append(img_band02)
        inImgBands.append(img_band03)
        inImgBands.append(img_band04)
        inImgBands.append(img_band05)
        inImgBands.append(img_band06)
        inImgBands.append(img_band07)
        inImgBands.append(img_band08)
        inImgBands.append(img_band08A)
        inImgBands.append(img_band11)
        inImgBands.append(img_band12)

        solarIrrVals = list()
        IrrVal = collections.namedtuple("SolarIrradiance", ["irradiance"])
        solarIrrVals.append(IrrVal(irradiance=self.esun_B1))
        solarIrrVals.append(IrrVal(irradiance=self.esun_B2))
        solarIrrVals.append(IrrVal(irradiance=self.esun_B3))
        solarIrrVals.append(IrrVal(irradiance=self.esun_B4))
        solarIrrVals.append(IrrVal(irradiance=self.esun_B5))
        solarIrrVals.append(IrrVal(irradiance=self.esun_B6))
        solarIrrVals.append(IrrVal(irradiance=self.esun_B7))
        solarIrrVals.append(IrrVal(irradiance=self.esun_B8))
        solarIrrVals.append(IrrVal(irradiance=self.esun_B11))
        solarIrrVals.append(IrrVal(irradiance=self.esun_B12))

        rsgislib.imagecalibration.toa_refl_to_radiance(
            inImgBands,
            outputReflImage,
            outFormat,
            rsgislib.TYPE_32FLOAT,
            self.quantificationVal,
            self.earthSunDist_U,
            self.solarZenith,
            solarIrrVals,
        )

        if (img_band02_tmp is not None) and os.path.exists(img_band02_tmp):
            os.remove(img_band02_tmp)

        if (img_band03_tmp is not None) and os.path.exists(img_band03_tmp):
            os.remove(img_band03_tmp)

        if (img_band04_tmp is not None) and os.path.exists(img_band04_tmp):
            os.remove(img_band04_tmp)

        if (img_band05_tmp is not None) and os.path.exists(img_band05_tmp):
            os.remove(img_band05_tmp)

        if (img_band06_tmp is not None) and os.path.exists(img_band06_tmp):
            os.remove(img_band06_tmp)

        if (img_band07_tmp is not None) and os.path.exists(img_band07_tmp):
            os.remove(img_band07_tmp)

        if (img_band08_tmp is not None) and os.path.exists(img_band08_tmp):
            os.remove(img_band08_tmp)

        if (img_band08A_tmp is not None) and os.path.exists(img_band08A_tmp):
            os.remove(img_band08A_tmp)

        if (img_band11_tmp is not None) and os.path.exists(img_band11_tmp):
            os.remove(img_band11_tmp)

        if (img_band12_tmp is not None) and os.path.exists(img_band12_tmp):
            os.remove(img_band12_tmp)

        return outputReflImage, outputThermalImage

    def convertThermalToBrightness(
        self, inputRadImage, outputPath, outputName, outFormat, scaleFactor
    ):
        raise ARCSIException("Not Implemented")

    def convertImageToTOARefl(
        self, inputRadImage, outputPath, outputName, outFormat, scaleFactor
    ):
        print("Converting to TOA")
        self.imgIntScaleFactor = scaleFactor

        outputImage = os.path.join(outputPath, outputName)
        solarIrrVals = list()
        IrrVal = collections.namedtuple("SolarIrradiance", ["irradiance"])
        solarIrrVals.append(IrrVal(irradiance=self.esun_B1))
        solarIrrVals.append(IrrVal(irradiance=self.esun_B2))
        solarIrrVals.append(IrrVal(irradiance=self.esun_B3))
        solarIrrVals.append(IrrVal(irradiance=self.esun_B4))
        solarIrrVals.append(IrrVal(irradiance=self.esun_B5))
        solarIrrVals.append(IrrVal(irradiance=self.esun_B6))
        solarIrrVals.append(IrrVal(irradiance=self.esun_B7))
        solarIrrVals.append(IrrVal(irradiance=self.esun_B8))
        solarIrrVals.append(IrrVal(irradiance=self.esun_B11))
        solarIrrVals.append(IrrVal(irradiance=self.esun_B12))
        rsgislib.imagecalibration.radiance_to_toa_refl(
            inputRadImage,
            outputImage,
            outFormat,
            rsgislib.TYPE_16UINT,
            scaleFactor,
            self.acquisitionTime.year,
            self.acquisitionTime.month,
            self.acquisitionTime.day,
            self.solarZenith,
            solarIrrVals,
        )
        return outputImage

    def generateCloudMask(
        self,
        inputReflImage,
        inputSatImage,
        inputThermalImage,
        inputViewAngleImg,
        inputValidImg,
        outputPath,
        outputCloudName,
        outputCloudProb,
        outFormat,
        tmpPath,
        scaleFactor,
        cloud_msk_methods=None,
    ):
        outputImage = os.path.join(outputPath, outputCloudName)
        outputProbImage = None
        tmpBaseName = os.path.splitext(outputCloudName)[0]
        tmpBaseDIR = os.path.join(tmpPath, tmpBaseName)

        tmpDIRExisted = True
        if not os.path.exists(tmpBaseDIR):
            os.makedirs(tmpBaseDIR)
            tmpDIRExisted = False

        #########################################################################################################
        # Create 13 band tmp images which has all image bands for input into cloud masking.
        sen2ImgB01_tmp = os.path.join(tmpBaseDIR, f"{tmpBaseName}_B01.kea")
        rsgislib.imageutils.resample_img_to_match(
            inputReflImage,
            self.sen2ImgB01,
            sen2ImgB01_tmp,
            "KEA",
            rsgislib.INTERP_NEAREST_NEIGHBOUR,
            rsgislib.TYPE_16UINT,
            multicore=False,
        )
        sen2ImgB09_tmp = os.path.join(tmpBaseDIR, f"{tmpBaseName}_B09.kea")
        rsgislib.imageutils.resample_img_to_match(
            inputReflImage,
            self.sen2ImgB09,
            sen2ImgB09_tmp,
            "KEA",
            rsgislib.INTERP_NEAREST_NEIGHBOUR,
            rsgislib.TYPE_16UINT,
            multicore=False,
        )
        sen2ImgB10_tmp = os.path.join(tmpBaseDIR, f"{tmpBaseName}_B10.kea")
        rsgislib.imageutils.resample_img_to_match(
            inputReflImage,
            self.sen2ImgB10,
            sen2ImgB10_tmp,
            "KEA",
            rsgislib.INTERP_NEAREST_NEIGHBOUR,
            rsgislib.TYPE_16UINT,
            multicore=False,
        )

        if abs(self.ratiometric_offs_B0) > 0:
            basename = os.path.splitext(os.path.basename(sen2ImgB01_tmp))[0]
            sen2ImgB01_tmp2 = os.path.join(
                tmpBaseDIR, f"{basename}_tmpoffapply.kea"
            )
            exp = f"b1 + {self.ratiometric_offs_B0}"
            rsgislib.imagecalc.image_math(
                sen2ImgB01_tmp, sen2ImgB01_tmp2, exp, "KEA", rsgislib.TYPE_32FLOAT
            )
            sen2ImgB01_tmp = sen2ImgB01_tmp2

        if abs(self.ratiometric_offs_B9) > 0:
            basename = os.path.splitext(os.path.basename(sen2ImgB09_tmp))[0]
            sen2ImgB09_tmp2 = os.path.join(
                tmpBaseDIR, f"{basename}_tmpoffapply.kea"
            )
            exp = f"b1 + {self.ratiometric_offs_B9}"
            rsgislib.imagecalc.image_math(
                sen2ImgB09_tmp, sen2ImgB09_tmp2, exp, "KEA", rsgislib.TYPE_32FLOAT
            )
            sen2ImgB09_tmp = sen2ImgB09_tmp2

        if abs(self.ratiometric_offs_B10) > 0:
            basename = os.path.splitext(os.path.basename(sen2ImgB10_tmp))[0]
            sen2ImgB10_tmp2 = os.path.join(
                tmpBaseDIR, f"{basename}_tmpoffapply.kea"
            )
            exp = f"b1 + {self.ratiometric_offs_B10}"
            rsgislib.imagecalc.image_math(
                sen2ImgB10_tmp, sen2ImgB10_tmp2, exp, "KEA", rsgislib.TYPE_32FLOAT
            )
            sen2ImgB10_tmp = sen2ImgB10_tmp2

        # Stack Image Bands
        tmpTOAImg = os.path.join(tmpBaseDIR, f"{tmpBaseName}_pyfmasktmpTOA.kea")
        if self.imgIntScaleFactor != 10000:
            vrtImgB1B9B10 = os.path.join(
                tmpBaseDIR, f"{tmpBaseName}_b01b09b10_60mBands.vrt"
            )
            rsgislib.imageutils.create_stack_images_vrt(
                [sen2ImgB01_tmp, sen2ImgB09_tmp, sen2ImgB10_tmp], vrtImgB1B9B10
            )

            sen2ImgB1B9B10Rescaled = os.path.join(
                tmpBaseDIR, f"{tmpBaseName}_B01B09B10Rescaled.kea"
            )
            rsgislib.imagecalc.image_math(
                vrtImgB1B9B10,
                sen2ImgB1B9B10Rescaled,
                f"(b1/10000)*{self.imgIntScaleFactor}",
                "KEA",
                rsgislib.TYPE_16UINT,
            )
            rsgislib.imageutils.stack_img_bands(
                [inputReflImage, sen2ImgB1B9B10Rescaled],
                None,
                tmpTOAImg,
                None,
                0,
                "KEA",
                rsgislib.TYPE_16UINT,
            )
        else:
            rsgislib.imageutils.stack_img_bands(
                [inputReflImage, sen2ImgB01_tmp, sen2ImgB09_tmp, sen2ImgB10_tmp],
                None,
                tmpTOAImg,
                None,
                0,
                "KEA",
                rsgislib.TYPE_16UINT,
            )

        # Re-order image bands to be in correct order.
        fmaskReflImg = os.path.join(tmpBaseDIR, f"{tmpBaseName}_pyfmaskRefl.kea")
        rsgislib.imageutils.select_img_bands(
            tmpTOAImg,
            fmaskReflImg,
            "KEA",
            rsgislib.TYPE_16UINT,
            [11, 1, 2, 3, 4, 5, 6, 7, 8, 12, 13, 9, 10],
        )
        #########################################################################################################

        if (cloud_msk_methods is None) or (
            (cloud_msk_methods == "FMASK") or (cloud_msk_methods == "FMASK_DISP")
        ):
            anglesInfo = fmask.config.AnglesFileInfo(
                inputViewAngleImg,
                3,
                inputViewAngleImg,
                2,
                inputViewAngleImg,
                1,
                inputViewAngleImg,
                0,
            )
            fmaskCloudsImg = os.path.join(
                tmpBaseDIR, f"{tmpBaseName}_pyfmaskCloudsResult.kea"
            )
            fmaskFilenames = fmask.config.FmaskFilenames()
            fmaskFilenames.setTOAReflectanceFile(fmaskReflImg)
            fmaskFilenames.setSaturationMask(inputSatImage)
            fmaskFilenames.setOutputCloudMaskFile(fmaskCloudsImg)

            fmaskConfig = fmask.config.FmaskConfig(fmask.config.FMASK_SENTINEL2)
            fmaskConfig.setAnglesInfo(anglesInfo)
            fmaskConfig.setKeepIntermediates(True)
            fmaskConfig.setVerbose(True)
            fmaskConfig.setTempDir(tmpBaseDIR)
            fmaskConfig.setTOARefScaling(float(self.imgIntScaleFactor))
            fmaskConfig.setMinCloudSize(8)
            off_dict = {
                fmask.config.BAND_BLUE: 0,
                fmask.config.BAND_GREEN: 0,
                fmask.config.BAND_RED: 0,
                fmask.config.BAND_NIR: 0,
                fmask.config.BAND_SWIR1: 0,
                fmask.config.BAND_SWIR2: 0,
                fmask.config.BAND_WATERVAPOUR: 0,
                fmask.config.BAND_CIRRUS: 0,
                fmask.config.BAND_S2CDI_NIR7: 0,
                fmask.config.BAND_S2CDI_NIR8A: 0,
            }
            fmaskConfig.setTOARefOffsetDict(off_dict)

            if (cloud_msk_methods is not None) and ("FMASK_DISP" in cloud_msk_methods):
                fmaskConfig.setSen2displacementTest(
                    True
                )  # Frantz et al implementation.
            else:
                fmaskConfig.setSen2displacementTest(
                    False
                )  # Frantz et al implementation.
            fmaskConfig.setCloudBufferSize(10)
            fmaskConfig.setShadowBufferSize(10)

            fmask.fmask.doFmask(fmaskFilenames, fmaskConfig)
            rsgislib.imagecalc.image_math(
                fmaskCloudsImg,
                outputImage,
                "(b1==2)?1:(b1==3)?2:0",
                outFormat,
                rsgislib.TYPE_8UINT,
            )
        elif "S2CLOUDLESS" in cloud_msk_methods:
            from arcsilib.s2cloudless import run_pyfmask_shadow_masking, run_s2cloudless

            outputProbImage = os.path.join(outputPath, outputCloudProb)
            out_cloud_msk = os.path.join(
                tmpBaseDIR, f"{tmpBaseName}_s2cloudless_cloud_msk.kea"
            )

            run_s2cloudless(
                fmaskReflImg,
                out_cloud_msk,
                outputProbImage,
                inputValidImg,
                outFormat,
                tmpBaseDIR,
                toa_scale_factor=float(self.imgIntScaleFactor),
                min_obj_size=10,
                morph_close_size=5,
                morph_dilate_size=9,
            )

            run_pyfmask_shadow_masking(
                fmaskReflImg,
                inputSatImage,
                inputViewAngleImg,
                out_cloud_msk,
                tmpBaseDIR,
                float(self.imgIntScaleFactor),
                outputImage,
            )
        elif ("S2LESSFMSK" in cloud_msk_methods) or (
            "S2LESSFMSKD" in cloud_msk_methods
        ):
            from arcsilib.s2cloudless import (
                run_fmask_cloud_msk,
                run_pyfmask_shadow_masking,
                run_s2cloudless,
            )

            out_s2less_cloud_msk = os.path.join(
                tmpBaseDIR, f"{tmpBaseName}_s2cloudless_cloud_msk.kea"
            )
            outputProbImage = os.path.join(outputPath, outputCloudProb)
            run_s2cloudless(
                fmaskReflImg,
                out_s2less_cloud_msk,
                outputProbImage,
                inputValidImg,
                outFormat,
                tmpBaseDIR,
                toa_scale_factor=float(self.imgIntScaleFactor),
                min_obj_size=10,
                morph_close_size=5,
                morph_dilate_size=9,
            )

            out_fmsk_cloud_msk = os.path.join(
                tmpBaseDIR, f"{tmpBaseName}_fmsk_cloud_msk.kea"
            )
            use_frantz_disp = False
            if "S2LESSFMSKD" in cloud_msk_methods:
                use_frantz_disp = True
            run_fmask_cloud_msk(
                fmaskReflImg,
                inputSatImage,
                inputViewAngleImg,
                out_fmsk_cloud_msk,
                tmpBaseDIR,
                float(self.imgIntScaleFactor),
                use_frantz_disp,
            )

            # Combine cloud masks
            out_cloud_msk = os.path.join(
                tmpBaseDIR, f"{tmpBaseName}_fmsk_s2l_cloud_msk.kea"
            )
            bandDefns = []
            bandDefns.append(
                rsgislib.imagecalc.BandDefn("s2l", out_s2less_cloud_msk, 1)
            )
            bandDefns.append(rsgislib.imagecalc.BandDefn("fmsk", out_fmsk_cloud_msk, 1))
            rsgislib.imagecalc.band_math(
                out_cloud_msk,
                "(s2l==1)&&(fmsk==1)?1:0",
                outFormat,
                rsgislib.TYPE_8UINT,
                bandDefns,
            )

            # Remove small cloud features
            out_cloud_msk_clumps = os.path.join(
                tmpBaseDIR, f"{tmpBaseName}_fmsk_s2l_cloud_msk_clumps.kea"
            )
            rsgislib.segmentation.clump(
                out_cloud_msk, out_cloud_msk_clumps, "KEA", False, 0, False
            )
            rsgislib.rastergis.pop_rat_img_stats(
                out_cloud_msk_clumps, True, False, True
            )
            out_cloud_msk_clumps_rmsml = os.path.join(
                tmpBaseDIR, f"{tmpBaseName}_fmsk_s2l_cloud_msk_clumps_rmsml.kea"
            )
            rsgislib.segmentation.rm_small_clumps(
                out_cloud_msk_clumps, out_cloud_msk_clumps_rmsml, 8, "KEA"
            )
            out_cloud_msk_rmsml = os.path.join(
                tmpBaseDIR, f"{tmpBaseName}_fmsk_s2l_cloud_msk_rmsml.kea"
            )
            rsgislib.imagecalc.image_math(
                out_cloud_msk_clumps_rmsml,
                out_cloud_msk_rmsml,
                "b1>0?1:0",
                "KEA",
                rsgislib.TYPE_8UINT,
            )

            # Buffer the cloud features
            morph_operator = os.path.join(tmpBaseDIR, "morph_circ5")
            morph_operator_file = f"{morph_operator}.gmtxt"
            morph_op_size = 5
            rsgislib.imagemorphology.create_circular_op(morph_operator, morph_op_size)
            out_cloud_msk_dil = os.path.join(
                tmpBaseDIR, f"{tmpBaseName}_fmsk_s2l_cloud_msk_dilate.kea"
            )
            rsgislib.imagemorphology.image_dilate(
                out_cloud_msk_rmsml,
                out_cloud_msk_dil,
                morph_operator_file,
                True,
                5,
                "KEA",
                rsgislib.TYPE_8UINT,
            )

            # Find shadow mask and add it to cloud mask.
            run_pyfmask_shadow_masking(
                fmaskReflImg,
                inputSatImage,
                inputViewAngleImg,
                out_cloud_msk_dil,
                tmpBaseDIR,
                float(self.imgIntScaleFactor),
                outputImage,
            )
        else:
            raise ARCSIException(
                "The cloud masking methods specified for Sentinel-2 does not exist."
            )

        if outFormat == "KEA":
            rsgislib.rastergis.pop_rat_img_stats(outputImage, True, True, True)
            ratDataset = gdal.Open(outputImage, gdal.GA_Update)
            red = rat.readColumn(ratDataset, "Red")
            green = rat.readColumn(ratDataset, "Green")
            blue = rat.readColumn(ratDataset, "Blue")
            ClassName = numpy.empty_like(red, dtype=numpy.dtype("a255"))

            red[0] = 0
            green[0] = 0
            blue[0] = 0

            if (red.shape[0] == 2) or (red.shape[0] == 3):
                red[1] = 0
                green[1] = 0
                blue[1] = 255
                ClassName[1] = "Clouds"

                if red.shape[0] == 3:
                    red[2] = 0
                    green[2] = 255
                    blue[2] = 255
                    ClassName[2] = "Shadows"

            rat.writeColumn(ratDataset, "Red", red)
            rat.writeColumn(ratDataset, "Green", green)
            rat.writeColumn(ratDataset, "Blue", blue)
            rat.writeColumn(ratDataset, "ClassName", ClassName)
            ratDataset = None
        rsgislib.imageutils.copy_proj_from_img(outputImage, inputReflImage)
        if outputProbImage is not None:
            rsgislib.imageutils.copy_proj_from_img(outputProbImage, inputReflImage)

        if not self.debugMode:
            if not tmpDIRExisted:
                shutil.rmtree(tmpBaseDIR, ignore_errors=True)

        return outputImage, outputProbImage

    def createCloudMaskDataArray(self, inImgDataArr):
        return inImgDataArr

    def defineDarkShadowImageBand(self):
        return 7

    def calc6SCoefficients(
        self, aeroProfile, atmosProfile, grdRefl, surfaceAltitude, aotVal, useBRDF
    ):
        sixsCoeffs = numpy.zeros((10, 6), dtype=numpy.float32)
        # Set up 6S model
        s = Py6S.SixS()
        s.atmos_profile = atmosProfile
        s.aero_profile = aeroProfile
        s.ground_reflectance = grdRefl
        s.geometry = Py6S.Geometry.User()
        s.geometry.solar_z = self.solarZenith
        s.geometry.solar_a = self.solarAzimuth
        s.geometry.view_z = self.sensorZenith
        s.geometry.view_a = self.sensorAzimuth
        s.geometry.month = self.acquisitionTime.month
        s.geometry.day = self.acquisitionTime.day
        s.geometry.gmt_decimal_hour = (
            float(self.acquisitionTime.hour) + float(self.acquisitionTime.minute) / 60.0
        )
        s.geometry.latitude = self.latCentre
        s.geometry.longitude = self.lonCentre
        s.altitudes = Py6S.Altitudes()
        s.altitudes.set_target_custom_altitude(surfaceAltitude)
        s.altitudes.set_sensor_satellite_level()
        if useBRDF:
            s.atmos_corr = Py6S.AtmosCorr.AtmosCorrBRDFFromRadiance(200)
        else:
            s.atmos_corr = Py6S.AtmosCorr.AtmosCorrLambertianFromRadiance(200)
        s.aot550 = aotVal

        # Blue
        s.wavelength = Py6S.Wavelength(
            self.specBandInfo["B2"].wvLenMin6S,
            self.specBandInfo["B2"].wvLenMax6S,
            self.specBandInfo["B2"].respFunc6S,
        )
        s.run()
        sixsCoeffs[0, 0] = float(s.outputs.values["coef_xa"])
        sixsCoeffs[0, 1] = float(s.outputs.values["coef_xb"])
        sixsCoeffs[0, 2] = float(s.outputs.values["coef_xc"])
        sixsCoeffs[0, 3] = float(s.outputs.values["direct_solar_irradiance"])
        sixsCoeffs[0, 4] = float(s.outputs.values["diffuse_solar_irradiance"])
        sixsCoeffs[0, 5] = float(s.outputs.values["environmental_irradiance"])

        # Green
        s.wavelength = Py6S.Wavelength(
            self.specBandInfo["B3"].wvLenMin6S,
            self.specBandInfo["B3"].wvLenMax6S,
            self.specBandInfo["B3"].respFunc6S,
        )
        s.run()
        sixsCoeffs[1, 0] = float(s.outputs.values["coef_xa"])
        sixsCoeffs[1, 1] = float(s.outputs.values["coef_xb"])
        sixsCoeffs[1, 2] = float(s.outputs.values["coef_xc"])
        sixsCoeffs[1, 3] = float(s.outputs.values["direct_solar_irradiance"])
        sixsCoeffs[1, 4] = float(s.outputs.values["diffuse_solar_irradiance"])
        sixsCoeffs[1, 5] = float(s.outputs.values["environmental_irradiance"])

        # Red
        s.wavelength = Py6S.Wavelength(
            self.specBandInfo["B4"].wvLenMin6S,
            self.specBandInfo["B4"].wvLenMax6S,
            self.specBandInfo["B4"].respFunc6S,
        )
        s.run()
        sixsCoeffs[2, 0] = float(s.outputs.values["coef_xa"])
        sixsCoeffs[2, 1] = float(s.outputs.values["coef_xb"])
        sixsCoeffs[2, 2] = float(s.outputs.values["coef_xc"])
        sixsCoeffs[2, 3] = float(s.outputs.values["direct_solar_irradiance"])
        sixsCoeffs[2, 4] = float(s.outputs.values["diffuse_solar_irradiance"])
        sixsCoeffs[2, 5] = float(s.outputs.values["environmental_irradiance"])

        # RE B5
        s.wavelength = Py6S.Wavelength(
            self.specBandInfo["B5"].wvLenMin6S,
            self.specBandInfo["B5"].wvLenMax6S,
            self.specBandInfo["B5"].respFunc6S,
        )
        s.run()
        sixsCoeffs[3, 0] = float(s.outputs.values["coef_xa"])
        sixsCoeffs[3, 1] = float(s.outputs.values["coef_xb"])
        sixsCoeffs[3, 2] = float(s.outputs.values["coef_xc"])
        sixsCoeffs[3, 3] = float(s.outputs.values["direct_solar_irradiance"])
        sixsCoeffs[3, 4] = float(s.outputs.values["diffuse_solar_irradiance"])
        sixsCoeffs[3, 5] = float(s.outputs.values["environmental_irradiance"])

        # RE B6
        s.wavelength = Py6S.Wavelength(
            self.specBandInfo["B6"].wvLenMin6S,
            self.specBandInfo["B6"].wvLenMax6S,
            self.specBandInfo["B6"].respFunc6S,
        )
        s.run()
        sixsCoeffs[4, 0] = float(s.outputs.values["coef_xa"])
        sixsCoeffs[4, 1] = float(s.outputs.values["coef_xb"])
        sixsCoeffs[4, 2] = float(s.outputs.values["coef_xc"])
        sixsCoeffs[4, 3] = float(s.outputs.values["direct_solar_irradiance"])
        sixsCoeffs[4, 4] = float(s.outputs.values["diffuse_solar_irradiance"])
        sixsCoeffs[4, 5] = float(s.outputs.values["environmental_irradiance"])

        # RE B7
        s.wavelength = Py6S.Wavelength(
            self.specBandInfo["B7"].wvLenMin6S,
            self.specBandInfo["B7"].wvLenMax6S,
            self.specBandInfo["B7"].respFunc6S,
        )
        s.run()
        sixsCoeffs[5, 0] = float(s.outputs.values["coef_xa"])
        sixsCoeffs[5, 1] = float(s.outputs.values["coef_xb"])
        sixsCoeffs[5, 2] = float(s.outputs.values["coef_xc"])
        sixsCoeffs[5, 3] = float(s.outputs.values["direct_solar_irradiance"])
        sixsCoeffs[5, 4] = float(s.outputs.values["diffuse_solar_irradiance"])
        sixsCoeffs[5, 5] = float(s.outputs.values["environmental_irradiance"])

        # NIR B8
        s.wavelength = Py6S.Wavelength(
            self.specBandInfo["B8"].wvLenMin6S,
            self.specBandInfo["B8"].wvLenMax6S,
            self.specBandInfo["B8"].respFunc6S,
        )
        s.run()
        sixsCoeffs[6, 0] = float(s.outputs.values["coef_xa"])
        sixsCoeffs[6, 1] = float(s.outputs.values["coef_xb"])
        sixsCoeffs[6, 2] = float(s.outputs.values["coef_xc"])
        sixsCoeffs[6, 3] = float(s.outputs.values["direct_solar_irradiance"])
        sixsCoeffs[6, 4] = float(s.outputs.values["diffuse_solar_irradiance"])
        sixsCoeffs[6, 5] = float(s.outputs.values["environmental_irradiance"])

        # NIR B8A
        s.wavelength = Py6S.Wavelength(
            self.specBandInfo["B8A"].wvLenMin6S,
            self.specBandInfo["B8A"].wvLenMax6S,
            self.specBandInfo["B8A"].respFunc6S,
        )
        s.run()
        sixsCoeffs[7, 0] = float(s.outputs.values["coef_xa"])
        sixsCoeffs[7, 1] = float(s.outputs.values["coef_xb"])
        sixsCoeffs[7, 2] = float(s.outputs.values["coef_xc"])
        sixsCoeffs[7, 3] = float(s.outputs.values["direct_solar_irradiance"])
        sixsCoeffs[7, 4] = float(s.outputs.values["diffuse_solar_irradiance"])
        sixsCoeffs[7, 5] = float(s.outputs.values["environmental_irradiance"])

        # SWIR 1
        s.wavelength = Py6S.Wavelength(
            self.specBandInfo["B11"].wvLenMin6S,
            self.specBandInfo["B11"].wvLenMax6S,
            self.specBandInfo["B11"].respFunc6S,
        )
        s.run()
        sixsCoeffs[8, 0] = float(s.outputs.values["coef_xa"])
        sixsCoeffs[8, 1] = float(s.outputs.values["coef_xb"])
        sixsCoeffs[8, 2] = float(s.outputs.values["coef_xc"])
        sixsCoeffs[8, 3] = float(s.outputs.values["direct_solar_irradiance"])
        sixsCoeffs[8, 4] = float(s.outputs.values["diffuse_solar_irradiance"])
        sixsCoeffs[8, 5] = float(s.outputs.values["environmental_irradiance"])

        # SWIR 2
        s.wavelength = Py6S.Wavelength(
            self.specBandInfo["B12"].wvLenMin6S,
            self.specBandInfo["B12"].wvLenMax6S,
            self.specBandInfo["B12"].respFunc6S,
        )
        s.run()
        sixsCoeffs[9, 0] = float(s.outputs.values["coef_xa"])
        sixsCoeffs[9, 1] = float(s.outputs.values["coef_xb"])
        sixsCoeffs[9, 2] = float(s.outputs.values["coef_xc"])
        sixsCoeffs[9, 3] = float(s.outputs.values["direct_solar_irradiance"])
        sixsCoeffs[9, 4] = float(s.outputs.values["diffuse_solar_irradiance"])
        sixsCoeffs[9, 5] = float(s.outputs.values["environmental_irradiance"])

        return sixsCoeffs

    def convertImageToSurfaceReflSglParam(
        self,
        inputRadImage,
        outputPath,
        outputName,
        outFormat,
        aeroProfile,
        atmosProfile,
        grdRefl,
        surfaceAltitude,
        aotVal,
        useBRDF,
        scaleFactor,
    ):
        print("Converting to Surface Reflectance")
        outputImage = os.path.join(outputPath, outputName)

        imgBandCoeffs = list()

        sixsCoeffs = self.calc6SCoefficients(
            aeroProfile, atmosProfile, grdRefl, surfaceAltitude, aotVal, useBRDF
        )

        imgBandCoeffs.append(
            rsgislib.imagecalibration.Band6SCoeff(
                band=1,
                aX=float(sixsCoeffs[0, 0]),
                bX=float(sixsCoeffs[0, 1]),
                cX=float(sixsCoeffs[0, 2]),
                DirIrr=float(sixsCoeffs[0, 3]),
                DifIrr=float(sixsCoeffs[0, 4]),
                EnvIrr=float(sixsCoeffs[0, 5]),
            )
        )
        imgBandCoeffs.append(
            rsgislib.imagecalibration.Band6SCoeff(
                band=2,
                aX=float(sixsCoeffs[1, 0]),
                bX=float(sixsCoeffs[1, 1]),
                cX=float(sixsCoeffs[1, 2]),
                DirIrr=float(sixsCoeffs[1, 3]),
                DifIrr=float(sixsCoeffs[1, 4]),
                EnvIrr=float(sixsCoeffs[1, 5]),
            )
        )
        imgBandCoeffs.append(
            rsgislib.imagecalibration.Band6SCoeff(
                band=3,
                aX=float(sixsCoeffs[2, 0]),
                bX=float(sixsCoeffs[2, 1]),
                cX=float(sixsCoeffs[2, 2]),
                DirIrr=float(sixsCoeffs[2, 3]),
                DifIrr=float(sixsCoeffs[2, 4]),
                EnvIrr=float(sixsCoeffs[2, 5]),
            )
        )
        imgBandCoeffs.append(
            rsgislib.imagecalibration.Band6SCoeff(
                band=4,
                aX=float(sixsCoeffs[3, 0]),
                bX=float(sixsCoeffs[3, 1]),
                cX=float(sixsCoeffs[3, 2]),
                DirIrr=float(sixsCoeffs[3, 3]),
                DifIrr=float(sixsCoeffs[3, 4]),
                EnvIrr=float(sixsCoeffs[3, 5]),
            )
        )
        imgBandCoeffs.append(
            rsgislib.imagecalibration.Band6SCoeff(
                band=5,
                aX=float(sixsCoeffs[4, 0]),
                bX=float(sixsCoeffs[4, 1]),
                cX=float(sixsCoeffs[4, 2]),
                DirIrr=float(sixsCoeffs[4, 3]),
                DifIrr=float(sixsCoeffs[4, 4]),
                EnvIrr=float(sixsCoeffs[4, 5]),
            )
        )
        imgBandCoeffs.append(
            rsgislib.imagecalibration.Band6SCoeff(
                band=6,
                aX=float(sixsCoeffs[5, 0]),
                bX=float(sixsCoeffs[5, 1]),
                cX=float(sixsCoeffs[5, 2]),
                DirIrr=float(sixsCoeffs[5, 3]),
                DifIrr=float(sixsCoeffs[5, 4]),
                EnvIrr=float(sixsCoeffs[5, 5]),
            )
        )
        imgBandCoeffs.append(
            rsgislib.imagecalibration.Band6SCoeff(
                band=7,
                aX=float(sixsCoeffs[6, 0]),
                bX=float(sixsCoeffs[6, 1]),
                cX=float(sixsCoeffs[6, 2]),
                DirIrr=float(sixsCoeffs[6, 3]),
                DifIrr=float(sixsCoeffs[6, 4]),
                EnvIrr=float(sixsCoeffs[6, 5]),
            )
        )
        imgBandCoeffs.append(
            rsgislib.imagecalibration.Band6SCoeff(
                band=8,
                aX=float(sixsCoeffs[7, 0]),
                bX=float(sixsCoeffs[7, 1]),
                cX=float(sixsCoeffs[7, 2]),
                DirIrr=float(sixsCoeffs[7, 3]),
                DifIrr=float(sixsCoeffs[7, 4]),
                EnvIrr=float(sixsCoeffs[7, 5]),
            )
        )
        imgBandCoeffs.append(
            rsgislib.imagecalibration.Band6SCoeff(
                band=9,
                aX=float(sixsCoeffs[8, 0]),
                bX=float(sixsCoeffs[8, 1]),
                cX=float(sixsCoeffs[8, 2]),
                DirIrr=float(sixsCoeffs[8, 3]),
                DifIrr=float(sixsCoeffs[8, 4]),
                EnvIrr=float(sixsCoeffs[8, 5]),
            )
        )
        imgBandCoeffs.append(
            rsgislib.imagecalibration.Band6SCoeff(
                band=10,
                aX=float(sixsCoeffs[9, 0]),
                bX=float(sixsCoeffs[9, 1]),
                cX=float(sixsCoeffs[9, 2]),
                DirIrr=float(sixsCoeffs[9, 3]),
                DifIrr=float(sixsCoeffs[9, 4]),
                EnvIrr=float(sixsCoeffs[9, 5]),
            )
        )

        rsgislib.imagecalibration.apply_6s_coeff_single_param(
            inputRadImage,
            outputImage,
            outFormat,
            rsgislib.TYPE_16UINT,
            scaleFactor,
            0,
            True,
            imgBandCoeffs,
        )
        return outputImage

    def convertImageToSurfaceReflDEMElevLUT(
        self,
        inputRadImage,
        inputDEMFile,
        outputPath,
        outputName,
        outFormat,
        aeroProfile,
        atmosProfile,
        grdRefl,
        aotVal,
        useBRDF,
        surfaceAltitudeMin,
        surfaceAltitudeMax,
        scaleFactor,
        elevCoeffs=None,
    ):
        print("Converting to Surface Reflectance")
        outputImage = os.path.join(outputPath, outputName)

        if elevCoeffs is None:
            print("Build an LUT for elevation values.")
            elev6SCoeffsLUT = self.buildElevation6SCoeffLUT(
                aeroProfile,
                atmosProfile,
                grdRefl,
                aotVal,
                useBRDF,
                surfaceAltitudeMin,
                surfaceAltitudeMax,
            )
            print("LUT has been built.")

            elevCoeffs = list()
            for elevLUT in elev6SCoeffsLUT:
                imgBandCoeffs = list()
                sixsCoeffs = elevLUT.Coeffs
                elevVal = elevLUT.Elev
                imgBandCoeffs.append(
                    rsgislib.imagecalibration.Band6SCoeff(
                        band=1,
                        aX=float(sixsCoeffs[0, 0]),
                        bX=float(sixsCoeffs[0, 1]),
                        cX=float(sixsCoeffs[0, 2]),
                        DirIrr=float(sixsCoeffs[0, 3]),
                        DifIrr=float(sixsCoeffs[0, 4]),
                        EnvIrr=float(sixsCoeffs[0, 5]),
                    )
                )
                imgBandCoeffs.append(
                    rsgislib.imagecalibration.Band6SCoeff(
                        band=2,
                        aX=float(sixsCoeffs[1, 0]),
                        bX=float(sixsCoeffs[1, 1]),
                        cX=float(sixsCoeffs[1, 2]),
                        DirIrr=float(sixsCoeffs[1, 3]),
                        DifIrr=float(sixsCoeffs[1, 4]),
                        EnvIrr=float(sixsCoeffs[1, 5]),
                    )
                )
                imgBandCoeffs.append(
                    rsgislib.imagecalibration.Band6SCoeff(
                        band=3,
                        aX=float(sixsCoeffs[2, 0]),
                        bX=float(sixsCoeffs[2, 1]),
                        cX=float(sixsCoeffs[2, 2]),
                        DirIrr=float(sixsCoeffs[2, 3]),
                        DifIrr=float(sixsCoeffs[2, 4]),
                        EnvIrr=float(sixsCoeffs[2, 5]),
                    )
                )
                imgBandCoeffs.append(
                    rsgislib.imagecalibration.Band6SCoeff(
                        band=4,
                        aX=float(sixsCoeffs[3, 0]),
                        bX=float(sixsCoeffs[3, 1]),
                        cX=float(sixsCoeffs[3, 2]),
                        DirIrr=float(sixsCoeffs[3, 3]),
                        DifIrr=float(sixsCoeffs[3, 4]),
                        EnvIrr=float(sixsCoeffs[3, 5]),
                    )
                )
                imgBandCoeffs.append(
                    rsgislib.imagecalibration.Band6SCoeff(
                        band=5,
                        aX=float(sixsCoeffs[4, 0]),
                        bX=float(sixsCoeffs[4, 1]),
                        cX=float(sixsCoeffs[4, 2]),
                        DirIrr=float(sixsCoeffs[4, 3]),
                        DifIrr=float(sixsCoeffs[4, 4]),
                        EnvIrr=float(sixsCoeffs[4, 5]),
                    )
                )
                imgBandCoeffs.append(
                    rsgislib.imagecalibration.Band6SCoeff(
                        band=6,
                        aX=float(sixsCoeffs[5, 0]),
                        bX=float(sixsCoeffs[5, 1]),
                        cX=float(sixsCoeffs[5, 2]),
                        DirIrr=float(sixsCoeffs[5, 3]),
                        DifIrr=float(sixsCoeffs[5, 4]),
                        EnvIrr=float(sixsCoeffs[5, 5]),
                    )
                )
                imgBandCoeffs.append(
                    rsgislib.imagecalibration.Band6SCoeff(
                        band=7,
                        aX=float(sixsCoeffs[6, 0]),
                        bX=float(sixsCoeffs[6, 1]),
                        cX=float(sixsCoeffs[6, 2]),
                        DirIrr=float(sixsCoeffs[6, 3]),
                        DifIrr=float(sixsCoeffs[6, 4]),
                        EnvIrr=float(sixsCoeffs[6, 5]),
                    )
                )
                imgBandCoeffs.append(
                    rsgislib.imagecalibration.Band6SCoeff(
                        band=8,
                        aX=float(sixsCoeffs[7, 0]),
                        bX=float(sixsCoeffs[7, 1]),
                        cX=float(sixsCoeffs[7, 2]),
                        DirIrr=float(sixsCoeffs[7, 3]),
                        DifIrr=float(sixsCoeffs[7, 4]),
                        EnvIrr=float(sixsCoeffs[7, 5]),
                    )
                )
                imgBandCoeffs.append(
                    rsgislib.imagecalibration.Band6SCoeff(
                        band=9,
                        aX=float(sixsCoeffs[8, 0]),
                        bX=float(sixsCoeffs[8, 1]),
                        cX=float(sixsCoeffs[8, 2]),
                        DirIrr=float(sixsCoeffs[8, 3]),
                        DifIrr=float(sixsCoeffs[8, 4]),
                        EnvIrr=float(sixsCoeffs[8, 5]),
                    )
                )
                imgBandCoeffs.append(
                    rsgislib.imagecalibration.Band6SCoeff(
                        band=10,
                        aX=float(sixsCoeffs[9, 0]),
                        bX=float(sixsCoeffs[9, 1]),
                        cX=float(sixsCoeffs[9, 2]),
                        DirIrr=float(sixsCoeffs[9, 3]),
                        DifIrr=float(sixsCoeffs[9, 4]),
                        EnvIrr=float(sixsCoeffs[9, 5]),
                    )
                )
                elevCoeffs.append(
                    rsgislib.imagecalibration.ElevLUTFeat(
                        Elev=float(elevVal), Coeffs=imgBandCoeffs
                    )
                )

        rsgislib.imagecalibration.apply_6s_coeff_elev_lut_param(
            inputRadImage,
            inputDEMFile,
            outputImage,
            outFormat,
            rsgislib.TYPE_16UINT,
            scaleFactor,
            0,
            True,
            elevCoeffs,
        )
        return outputImage, elevCoeffs

    def convertImageToSurfaceReflAOTDEMElevLUT(
        self,
        inputRadImage,
        inputDEMFile,
        inputAOTImage,
        outputPath,
        outputName,
        outFormat,
        aeroProfile,
        atmosProfile,
        grdRefl,
        useBRDF,
        surfaceAltitudeMin,
        surfaceAltitudeMax,
        aotMin,
        aotMax,
        scaleFactor,
        elevAOTCoeffs=None,
    ):
        print("Converting to Surface Reflectance")
        outputImage = os.path.join(outputPath, outputName)

        if elevAOTCoeffs is None:
            print("Build an LUT for elevation and AOT values.")
            elevAOT6SCoeffsLUT = self.buildElevationAOT6SCoeffLUT(
                aeroProfile,
                atmosProfile,
                grdRefl,
                useBRDF,
                surfaceAltitudeMin,
                surfaceAltitudeMax,
                aotMin,
                aotMax,
            )

            elevAOTCoeffs = list()
            for elevLUT in elevAOT6SCoeffsLUT:
                elevVal = elevLUT.Elev
                aotLUT = elevLUT.Coeffs
                aot6SCoeffsOut = list()
                for aotFeat in aotLUT:
                    sixsCoeffs = aotFeat.Coeffs
                    aotVal = aotFeat.AOT
                    imgBandCoeffs = list()
                    imgBandCoeffs.append(
                        rsgislib.imagecalibration.Band6SCoeff(
                            band=1,
                            aX=float(sixsCoeffs[0, 0]),
                            bX=float(sixsCoeffs[0, 1]),
                            cX=float(sixsCoeffs[0, 2]),
                            DirIrr=float(sixsCoeffs[0, 3]),
                            DifIrr=float(sixsCoeffs[0, 4]),
                            EnvIrr=float(sixsCoeffs[0, 5]),
                        )
                    )
                    imgBandCoeffs.append(
                        rsgislib.imagecalibration.Band6SCoeff(
                            band=2,
                            aX=float(sixsCoeffs[1, 0]),
                            bX=float(sixsCoeffs[1, 1]),
                            cX=float(sixsCoeffs[1, 2]),
                            DirIrr=float(sixsCoeffs[1, 3]),
                            DifIrr=float(sixsCoeffs[1, 4]),
                            EnvIrr=float(sixsCoeffs[1, 5]),
                        )
                    )
                    imgBandCoeffs.append(
                        rsgislib.imagecalibration.Band6SCoeff(
                            band=3,
                            aX=float(sixsCoeffs[2, 0]),
                            bX=float(sixsCoeffs[2, 1]),
                            cX=float(sixsCoeffs[2, 2]),
                            DirIrr=float(sixsCoeffs[2, 3]),
                            DifIrr=float(sixsCoeffs[2, 4]),
                            EnvIrr=float(sixsCoeffs[2, 5]),
                        )
                    )
                    imgBandCoeffs.append(
                        rsgislib.imagecalibration.Band6SCoeff(
                            band=4,
                            aX=float(sixsCoeffs[3, 0]),
                            bX=float(sixsCoeffs[3, 1]),
                            cX=float(sixsCoeffs[3, 2]),
                            DirIrr=float(sixsCoeffs[3, 3]),
                            DifIrr=float(sixsCoeffs[3, 4]),
                            EnvIrr=float(sixsCoeffs[3, 5]),
                        )
                    )
                    imgBandCoeffs.append(
                        rsgislib.imagecalibration.Band6SCoeff(
                            band=5,
                            aX=float(sixsCoeffs[4, 0]),
                            bX=float(sixsCoeffs[4, 1]),
                            cX=float(sixsCoeffs[4, 2]),
                            DirIrr=float(sixsCoeffs[4, 3]),
                            DifIrr=float(sixsCoeffs[4, 4]),
                            EnvIrr=float(sixsCoeffs[4, 5]),
                        )
                    )
                    imgBandCoeffs.append(
                        rsgislib.imagecalibration.Band6SCoeff(
                            band=6,
                            aX=float(sixsCoeffs[5, 0]),
                            bX=float(sixsCoeffs[5, 1]),
                            cX=float(sixsCoeffs[5, 2]),
                            DirIrr=float(sixsCoeffs[5, 3]),
                            DifIrr=float(sixsCoeffs[5, 4]),
                            EnvIrr=float(sixsCoeffs[5, 5]),
                        )
                    )
                    imgBandCoeffs.append(
                        rsgislib.imagecalibration.Band6SCoeff(
                            band=7,
                            aX=float(sixsCoeffs[6, 0]),
                            bX=float(sixsCoeffs[6, 1]),
                            cX=float(sixsCoeffs[6, 2]),
                            DirIrr=float(sixsCoeffs[6, 3]),
                            DifIrr=float(sixsCoeffs[6, 4]),
                            EnvIrr=float(sixsCoeffs[6, 5]),
                        )
                    )
                    imgBandCoeffs.append(
                        rsgislib.imagecalibration.Band6SCoeff(
                            band=8,
                            aX=float(sixsCoeffs[7, 0]),
                            bX=float(sixsCoeffs[7, 1]),
                            cX=float(sixsCoeffs[7, 2]),
                            DirIrr=float(sixsCoeffs[7, 3]),
                            DifIrr=float(sixsCoeffs[7, 4]),
                            EnvIrr=float(sixsCoeffs[7, 5]),
                        )
                    )
                    imgBandCoeffs.append(
                        rsgislib.imagecalibration.Band6SCoeff(
                            band=9,
                            aX=float(sixsCoeffs[8, 0]),
                            bX=float(sixsCoeffs[8, 1]),
                            cX=float(sixsCoeffs[8, 2]),
                            DirIrr=float(sixsCoeffs[8, 3]),
                            DifIrr=float(sixsCoeffs[8, 4]),
                            EnvIrr=float(sixsCoeffs[8, 5]),
                        )
                    )
                    imgBandCoeffs.append(
                        rsgislib.imagecalibration.Band6SCoeff(
                            band=10,
                            aX=float(sixsCoeffs[9, 0]),
                            bX=float(sixsCoeffs[9, 1]),
                            cX=float(sixsCoeffs[9, 2]),
                            DirIrr=float(sixsCoeffs[9, 3]),
                            DifIrr=float(sixsCoeffs[9, 4]),
                            EnvIrr=float(sixsCoeffs[9, 5]),
                        )
                    )
                    aot6SCoeffsOut.append(
                        rsgislib.imagecalibration.AOTLUTFeat(
                            AOT=float(aotVal), Coeffs=imgBandCoeffs
                        )
                    )
                elevAOTCoeffs.append(
                    rsgislib.imagecalibration.ElevLUTFeat(
                        Elev=float(elevVal), Coeffs=aot6SCoeffsOut
                    )
                )

        rsgislib.imagecalibration.apply_6s_coeff_elev_aot_lut_param(
            inputRadImage,
            inputDEMFile,
            inputAOTImage,
            outputImage,
            outFormat,
            rsgislib.TYPE_16UINT,
            scaleFactor,
            0,
            True,
            elevAOTCoeffs,
        )

    def run6SToOptimiseAODValue(
        self,
        aotVal,
        radBlueVal,
        predBlueVal,
        aeroProfile,
        atmosProfile,
        grdRefl,
        surfaceAltitude,
    ):
        """Used as part of the optimastion for identifying values of AOD"""
        print(
            "Testing AOD Val: ",
            aotVal,
        )
        s = Py6S.SixS()
        s.atmos_profile = atmosProfile
        s.aero_profile = aeroProfile
        s.ground_reflectance = grdRefl
        s.geometry = Py6S.Geometry.User()
        s.geometry.solar_z = self.solarZenith
        s.geometry.solar_a = self.solarAzimuth
        s.geometry.view_z = self.sensorZenith
        s.geometry.view_a = self.sensorAzimuth
        s.geometry.month = self.acquisitionTime.month
        s.geometry.day = self.acquisitionTime.day
        s.geometry.gmt_decimal_hour = (
            float(self.acquisitionTime.hour) + float(self.acquisitionTime.minute) / 60.0
        )
        s.geometry.latitude = self.latCentre
        s.geometry.longitude = self.lonCentre
        s.altitudes = Py6S.Altitudes()
        s.altitudes.set_target_custom_altitude(surfaceAltitude)
        s.altitudes.set_sensor_satellite_level()
        s.atmos_corr = Py6S.AtmosCorr.AtmosCorrLambertianFromRadiance(200)
        s.aot550 = aotVal

        # Band 1 (Blue!)
        s.wavelength = Py6S.Wavelength(
            self.specBandInfo["B2"].wvLenMin6S,
            self.specBandInfo["B2"].wvLenMax6S,
            self.specBandInfo["B2"].respFunc6S,
        )
        s.run()
        aX = float(s.outputs.values["coef_xa"])
        bX = float(s.outputs.values["coef_xb"])
        cX = float(s.outputs.values["coef_xc"])

        tmpVal = (aX * radBlueVal) - bX
        reflBlueVal = tmpVal / (1.0 + cX * tmpVal)
        outDist = math.sqrt(math.pow((reflBlueVal - predBlueVal), 2))
        print("\taX: ", aX, " bX: ", bX, " cX: ", cX, "     Dist = ", outDist)
        return outDist

    def findDDVTargets(self, inputTOAImage, outputPath, outputName, outFormat, tmpPath):
        raise ARCSIException("Not Implemented")

    def estimateImageToAODUsingDDV(
        self,
        inputRADImage,
        inputTOAImage,
        inputDEMFile,
        shadowMask,
        outputPath,
        outputName,
        outFormat,
        tmpPath,
        aeroProfile,
        atmosProfile,
        grdRefl,
        aotValMin,
        aotValMax,
    ):
        raise ARCSIException("Not Implemented")

    def estimateImageToAODUsingDOS(
        self,
        inputRADImage,
        inputTOAImage,
        inputDEMFile,
        shadowMask,
        outputPath,
        outputName,
        outFormat,
        tmpPath,
        aeroProfile,
        atmosProfile,
        grdRefl,
        aotValMin,
        aotValMax,
        globalDOS,
        simpleDOS,
        dosOutRefl,
    ):
        raise ARCSIException("Not Implemented")

    def estimateSingleAOTFromDOS(
        self,
        radianceImage,
        toaImage,
        inputDEMFile,
        tmpPath,
        outputName,
        outFormat,
        aeroProfile,
        atmosProfile,
        grdRefl,
        minAOT,
        maxAOT,
        dosOutRefl,
    ):
        try:
            return self.estimateSingleAOTFromDOSBandImpl(
                radianceImage,
                toaImage,
                inputDEMFile,
                tmpPath,
                outputName,
                outFormat,
                aeroProfile,
                atmosProfile,
                grdRefl,
                minAOT,
                maxAOT,
                dosOutRefl,
                1,
            )
        except Exception as e:
            raise

    def setBandNames(self, imageFile):
        dataset = gdal.Open(imageFile, gdal.GA_Update)
        if not dataset is None:
            dataset.GetRasterBand(1).SetDescription("Blue")
            dataset.GetRasterBand(2).SetDescription("Green")
            dataset.GetRasterBand(3).SetDescription("Red")
            dataset.GetRasterBand(4).SetDescription("RE_B5")
            dataset.GetRasterBand(5).SetDescription("RE_B6")
            dataset.GetRasterBand(6).SetDescription("RE_B7")
            dataset.GetRasterBand(7).SetDescription("NIR_B8")
            dataset.GetRasterBand(8).SetDescription("NIR_B8A")
            dataset.GetRasterBand(9).SetDescription("SWIR1")
            dataset.GetRasterBand(10).SetDescription("SWIR2")
            dataset = None
        else:
            print("Could not open image to set band names: ", imageFile)

    def cleanLocalFollowProcessing(self):
        if not self.debugMode:
            gdalDriver = gdal.GetDriverByName("KEA")
            if self.resampleTo20m:
                gdalDriver.Delete(self.sen2ImgB02_20m)
                gdalDriver.Delete(self.sen2ImgB03_20m)
                gdalDriver.Delete(self.sen2ImgB04_20m)
                gdalDriver.Delete(self.sen2ImgB08_20m)
            else:
                gdalDriver.Delete(self.sen2ImgB05_10m)
                gdalDriver.Delete(self.sen2ImgB06_10m)
                gdalDriver.Delete(self.sen2ImgB07_10m)
                gdalDriver.Delete(self.sen2ImgB8A_10m)
                gdalDriver.Delete(self.sen2ImgB11_10m)
                gdalDriver.Delete(self.sen2ImgB12_10m)
