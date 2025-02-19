"""
Module that contains the ARCSILandsatTMSensor class.
"""
############################################################################
#  arcsisensorlandsat.py
#
#  Copyright 2013 ARCSI.
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
# Purpose:  A class for read the landsat sensor header file and applying
#           the pre-processing operations within ARCSI to the landsat 5 TM
#           datasets.
#
# Author: Pete Bunting
# Email: pfb@aber.ac.uk
# Date: 05/07/2013
# Version: 1.0
#
# History:
# Version 1.0 - Created.
#
############################################################################

import collections
import datetime
import json
import math
import os
import shutil

import fmask.config
import fmask.fmask
import fmask.landsatangles
import numpy
import Py6S
import rios.fileinfo
import rsgislib
import rsgislib.imagecalc
import rsgislib.imagecalibration
import rsgislib.imagecalibration.solarangles
import rsgislib.imageutils
import rsgislib.rastergis
import rsgislib.segmentation
import rsgislib.segmentation.shepherdseg
import rsgislib.tools.geometrytools
import rsgislib.tools.utils
from osgeo import gdal, osr
from rios import rat

from .arcsiexception import ARCSIException
from .arcsisensor import ARCSIAbstractSensor
from .arcsiutils import ARCSILandsatMetaUtils


class ARCSILandsatTMSensor(ARCSIAbstractSensor):
    """
    A class which represents the landsat 5 TM sensor to read
    header parameters and apply data processing operations.
    """

    def __init__(self, debugMode, inputImage):
        ARCSIAbstractSensor.__init__(self, debugMode, inputImage)
        self.sensor = "LS_TM"
        self.collection_num = 0

        self.band1File = ""
        self.band2File = ""
        self.band3File = ""
        self.band4File = ""
        self.band5File = ""
        self.band6File = ""
        self.band7File = ""
        self.bandQAFile = ""
        self.row = 0
        self.path = 0

        self.b1CalMin = 0
        self.b1CalMax = 0
        self.b2CalMin = 0
        self.b2CalMax = 0
        self.b3CalMin = 0
        self.b3CalMax = 0
        self.b4CalMin = 0
        self.b4CalMax = 0
        self.b5CalMin = 0
        self.b5CalMax = 0
        self.b6CalMin = 0
        self.b6CalMax = 0
        self.b7CalMin = 0
        self.b7CalMax = 0

        self.b1MinRad = 0.0
        self.b1MaxRad = 0.0
        self.b2MinRad = 0.0
        self.b2MaxRad = 0.0
        self.b3MinRad = 0.0
        self.b3MaxRad = 0.0
        self.b4MinRad = 0.0
        self.b4MaxRad = 0.0
        self.b5MinRad = 0.0
        self.b5MaxRad = 0.0
        self.b6MinRad = 0.0
        self.b6MaxRad = 0.0
        self.b7MinRad = 0.0
        self.b7MaxRad = 0.0

        self.sensorID = ""
        self.spacecraftID = ""
        self.cloudCover = 0.0
        self.cloudCoverLand = 0.0
        self.earthSunDistance = 0.0
        self.gridCellSizeRefl = 0.0
        self.gridCellSizeTherm = 0.0

    def extractHeaderParameters(self, inputHeader, wktStr):
        """
        Understands and parses the Landsat MTL header files
        """
        try:
            if not self.userSpInputImage is None:
                raise ARCSIException(
                    "Landsat sensor cannot accept a user specified image file - only the images in the header file will be used."
                )
            self.headerFileName = os.path.split(inputHeader)[1]

            print("Reading header file")
            hFile = open(inputHeader, "r")
            headerParams = dict()
            for line in hFile:
                line = line.strip()
                if line:
                    lineVals = line.split("=")
                    if len(lineVals) == 2:
                        if (lineVals[0].strip() != "GROUP") or (
                            lineVals[0].strip() != "END_GROUP"
                        ):
                            headerParams[lineVals[0].strip()] = (
                                lineVals[1].strip().replace('"', "")
                            )
            hFile.close()
            print("Extracting Header Values")
            # Get the sensor info.
            if (
                (headerParams["SPACECRAFT_ID"].upper() == "LANDSAT_5")
                or (headerParams["SPACECRAFT_ID"].upper() == "LANDSAT5")
            ) and (headerParams["SENSOR_ID"].upper() == "TM"):
                self.sensor = "LS5TM"
            elif (
                (headerParams["SPACECRAFT_ID"].upper() == "LANDSAT_4")
                or (headerParams["SPACECRAFT_ID"].upper() == "LANDSAT4")
            ) and (headerParams["SENSOR_ID"].upper() == "TM"):
                self.sensor = "LS4TM"
            else:
                raise ARCSIException(
                    "Do no recognise the spacecraft and sensor or combination."
                )

            self.sensorID = headerParams["SENSOR_ID"]
            self.spacecraftID = headerParams["SPACECRAFT_ID"]

            if headerParams["COLLECTION_NUMBER"] == "01":
                self.collection_num = 1
            elif headerParams["COLLECTION_NUMBER"] == "02":
                self.collection_num = 2
            else:
                raise ARCSIException(
                    "Can only process collection 1 and 2 data: {}".format(
                        headerParams["COLLECTION_NUMBER"]
                    )
                )

            # Get row/path
            try:
                self.row = int(headerParams["WRS_ROW"])
            except KeyError:
                self.row = int(headerParams["STARTING_ROW"])
            self.path = int(headerParams["WRS_PATH"])

            # Get date and time of the acquisition
            try:
                acData = headerParams["DATE_ACQUIRED"].split("-")
            except KeyError:
                acData = headerParams["ACQUISITION_DATE"].split("-")
            try:
                acTime = headerParams["SCENE_CENTER_TIME"].split(":")
            except KeyError:
                acTime = headerParams["SCENE_CENTER_SCAN_TIME"].split(":")

            secsTime = acTime[2].split(".")
            self.acquisitionTime = datetime.datetime(
                int(acData[0]),
                int(acData[1]),
                int(acData[2]),
                int(acTime[0]),
                int(acTime[1]),
                int(secsTime[0]),
            )

            self.solarZenith = 90 - rsgislib.tools.utils.str_to_float(
                headerParams["SUN_ELEVATION"]
            )
            self.solarAzimuth = rsgislib.tools.utils.str_to_float(
                headerParams["SUN_AZIMUTH"]
            )

            # Get the geographic lat/long corners of the image.
            geoCorners = ARCSILandsatMetaUtils.getGeographicCorners(headerParams)

            self.latTL = geoCorners[0]
            self.lonTL = geoCorners[1]
            self.latTR = geoCorners[2]
            self.lonTR = geoCorners[3]
            self.latBL = geoCorners[4]
            self.lonBL = geoCorners[5]
            self.latBR = geoCorners[6]
            self.lonBR = geoCorners[7]

            # Get the projected X/Y corners of the image
            projectedCorners = ARCSILandsatMetaUtils.getProjectedCorners(headerParams)

            self.xTL = projectedCorners[0]
            self.yTL = projectedCorners[1]
            self.xTR = projectedCorners[2]
            self.yTR = projectedCorners[3]
            self.xBL = projectedCorners[4]
            self.yBL = projectedCorners[5]
            self.xBR = projectedCorners[6]
            self.yBR = projectedCorners[7]

            # Get projection
            inProj = osr.SpatialReference()
            if headerParams["MAP_PROJECTION"] == "UTM":
                try:
                    datum = headerParams["DATUM"]
                    if datum != "WGS84":
                        raise ARCSIException(
                            "Datum not recogised. Expected 'WGS84' got '{}'".format(
                                datum
                            )
                        )
                except KeyError:
                    pass
                try:
                    ellipsoid = headerParams["ELLIPSOID"]
                    if ellipsoid != "WGS84":
                        raise ARCSIException(
                            "Ellipsoid not recogised. Expected 'WGS84' got '{}'".format(
                                ellipsoid
                            )
                        )
                except KeyError:
                    pass
                try:
                    utmZone = int(headerParams["UTM_ZONE"])
                except KeyError:
                    utmZone = int(headerParams["ZONE_NUMBER"])
                # FIXME: should this be hardcoded to north?
                utmCode = "WGS84UTM" + str(utmZone) + str("N")
                inProj.ImportFromEPSG(self.epsgCodes[utmCode])
            elif (
                (headerParams["MAP_PROJECTION"] == "PS")
                and (headerParams["DATUM"] == "WGS84")
                and (headerParams["ELLIPSOID"] == "WGS84")
            ):
                inProj.ImportFromWkt(
                    'PROJCS["PS WGS84", GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563, AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433],AUTHORITY["EPSG","4326"]],PROJECTION["Polar_Stereographic"],PARAMETER["latitude_of_origin",-71],PARAMETER["central_meridian",0],PARAMETER["scale_factor",1],PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["metre",1,AUTHORITY["EPSG","9001"]]]'
                )
            else:
                raise ARCSIException(
                    "Expecting Landsat to be projected in UTM or PolarStereographic (PS) with datum=WGS84 and ellipsoid=WGS84."
                )

            if self.inWKT == "":
                self.inWKT = inProj.ExportToWkt()

            # Check image is square!
            if not (
                (self.xTL == self.xBL)
                and (self.yTL == self.yTR)
                and (self.xTR == self.xBR)
                and (self.yBL == self.yBR)
            ):
                raise ARCSIException("Image is not square in projected coordinates.")

            self.xCentre = self.xTL + ((self.xTR - self.xTL) / 2)
            self.yCentre = self.yBR + ((self.yTL - self.yBR) / 2)

            (
                self.lonCentre,
                self.latCentre,
            ) = rsgislib.tools.geometrytools.reproj_point_to_wgs84(
                inProj, self.xCentre, self.yCentre
            )

            # print("Lat: " + str(self.latCentre) + " Long: " + str(self.lonCentre))

            metaFilenames = ARCSILandsatMetaUtils.getBandFilenames(headerParams, 7)

            filesDIR = os.path.dirname(inputHeader)

            self.band1File = os.path.join(filesDIR, metaFilenames[0])
            self.band2File = os.path.join(filesDIR, metaFilenames[1])
            self.band3File = os.path.join(filesDIR, metaFilenames[2])
            self.band4File = os.path.join(filesDIR, metaFilenames[3])
            self.band5File = os.path.join(filesDIR, metaFilenames[4])
            self.band6File = os.path.join(filesDIR, metaFilenames[5])
            self.band7File = os.path.join(filesDIR, metaFilenames[6])

            if "FILE_NAME_BAND_QUALITY" in headerParams:
                self.bandQAFile = os.path.join(
                    filesDIR, headerParams["FILE_NAME_BAND_QUALITY"]
                )
            elif "FILE_NAME_QUALITY_L1_PIXEL" in headerParams:
                self.bandQAFile = os.path.join(
                    filesDIR, headerParams["FILE_NAME_QUALITY_L1_PIXEL"]
                )
            else:
                print(
                    "Warning - the quality band is not available. Are you using collection 1 or 2 data?"
                )
                self.bandQAFile = ""

            metaQCalMinList = []
            metaQCalMaxList = []

            for i in range(1, 8):
                try:
                    metaQCalMinList.append(
                        rsgislib.tools.utils.str_to_float(
                            headerParams["QUANTIZE_CAL_MIN_BAND_{}".format(i)], 1.0
                        )
                    )
                    metaQCalMaxList.append(
                        rsgislib.tools.utils.str_to_float(
                            headerParams["QUANTIZE_CAL_MAX_BAND_{}".format(i)], 255.0
                        )
                    )
                except KeyError:
                    metaQCalMinList.append(
                        rsgislib.tools.utils.str_to_float(
                            headerParams["QCALMIN_BAND{}".format(i)], 1.0
                        )
                    )
                    metaQCalMaxList.append(
                        rsgislib.tools.utils.str_to_float(
                            headerParams["QCALMAX_BAND{}".format(i)], 255.0
                        )
                    )

            self.b1CalMin = metaQCalMinList[0]
            self.b1CalMax = metaQCalMaxList[0]
            self.b2CalMin = metaQCalMinList[1]
            self.b2CalMax = metaQCalMaxList[1]
            self.b3CalMin = metaQCalMinList[2]
            self.b3CalMax = metaQCalMaxList[2]
            self.b4CalMin = metaQCalMinList[3]
            self.b4CalMax = metaQCalMaxList[3]
            self.b5CalMin = metaQCalMinList[4]
            self.b5CalMax = metaQCalMaxList[4]
            self.b6CalMin = metaQCalMinList[5]
            self.b6CalMax = metaQCalMaxList[5]
            self.b7CalMin = metaQCalMinList[6]
            self.b7CalMax = metaQCalMaxList[6]

            lMin = [-1.520, -2.840, -1.170, -1.510, -0.370, 1.238, -0.150]
            lMax = [193.000, 365.000, 264.000, 221.000, 30.200, 15.303, 16.500]
            metaRadMinList = []
            metaRadMaxList = []
            for i in range(1, 8):
                try:
                    metaRadMinList.append(
                        rsgislib.tools.utils.str_to_float(
                            headerParams["RADIANCE_MINIMUM_BAND_{}".format(i)],
                            lMin[i - 1],
                        )
                    )
                    metaRadMaxList.append(
                        rsgislib.tools.utils.str_to_float(
                            headerParams["RADIANCE_MAXIMUM_BAND_{}".format(i)],
                            lMax[i - 1],
                        )
                    )
                except KeyError:
                    metaRadMinList.append(
                        rsgislib.tools.utils.str_to_float(
                            headerParams["LMIN_BAND{}".format(i)], lMin[i - 1]
                        )
                    )
                    metaRadMaxList.append(
                        rsgislib.tools.utils.str_to_float(
                            headerParams["LMAX_BAND{}".format(i)], lMax[i - 1]
                        )
                    )

            self.b1MinRad = metaRadMinList[0]
            self.b1MaxRad = metaRadMaxList[0]
            self.b2MinRad = metaRadMinList[1]
            self.b2MaxRad = metaRadMaxList[1]
            self.b3MinRad = metaRadMinList[2]
            self.b3MaxRad = metaRadMaxList[2]
            self.b4MinRad = metaRadMinList[3]
            self.b4MaxRad = metaRadMaxList[3]
            self.b5MinRad = metaRadMinList[4]
            self.b5MaxRad = metaRadMaxList[4]
            self.b6MinRad = metaRadMinList[5]
            self.b6MaxRad = metaRadMaxList[5]
            self.b7MinRad = metaRadMinList[6]
            self.b7MaxRad = metaRadMaxList[6]

            if "CLOUD_COVER" in headerParams:
                self.cloudCover = rsgislib.tools.utils.str_to_float(
                    headerParams["CLOUD_COVER"], 0.0
                )
            if "CLOUD_COVER_LAND" in headerParams:
                self.cloudCoverLand = rsgislib.tools.utils.str_to_float(
                    headerParams["CLOUD_COVER_LAND"], 0.0
                )
            if "EARTH_SUN_DISTANCE" in headerParams:
                self.earthSunDistance = rsgislib.tools.utils.str_to_float(
                    headerParams["EARTH_SUN_DISTANCE"], 0.0
                )
            if "GRID_CELL_SIZE_REFLECTIVE" in headerParams:
                self.gridCellSizeRefl = rsgislib.tools.utils.str_to_float(
                    headerParams["GRID_CELL_SIZE_REFLECTIVE"], 60.0
                )
            if "GRID_CELL_SIZE_THERMAL" in headerParams:
                self.gridCellSizeTherm = rsgislib.tools.utils.str_to_float(
                    headerParams["GRID_CELL_SIZE_THERMAL"], 30.0
                )

            if "FILE_DATE" in headerParams:
                fileDateStr = headerParams["FILE_DATE"].strip()
            else:
                fileDateStr = headerParams["DATE_PRODUCT_GENERATED"].strip()
            fileDateStr = fileDateStr.replace("Z", "")
            self.fileDateObj = datetime.datetime.strptime(
                fileDateStr, "%Y-%m-%dT%H:%M:%S"
            )

            # Read MTL header into python dict for python-fmask
            self.fmaskMTLInfo = fmask.config.readMTLFile(inputHeader)

        except Exception as e:
            raise e

    def getSolarIrrStdSolarGeom(self):
        """
        Get Solar Azimuth and Zenith as standard geometry.
        Azimuth: N=0, E=90, S=180, W=270.
        """
        solarAz = rsgislib.imagecalibration.solarangles.get_solar_irr_convention_solar_azimuth_from_usgs(
            self.solarAzimuth
        )
        return (solarAz, self.solarZenith)

    def getSensorViewGeom(self):
        """
        Get sensor viewing angles
        returns (viewAzimuth, viewZenith)
        """
        return (0.0, 0.0)

    def generateOutputBaseName(self):
        """
        Provides an implementation for the landsat sensor
        """
        rowpath = "r" + str(self.row) + "p" + str(self.path)
        outname = self.defaultGenBaseOutFileName()
        outname = outname + str("_") + rowpath
        return outname

    def generateMetaDataFile(
        self,
        outputPath,
        outputFileName,
        productsStr,
        validMaskImage="",
        footprintCalc=False,
        calcdValuesDict=None,
        outFilesDict=None,
    ):
        """
        Generate file metadata.
        """
        if outFilesDict is None:
            outFilesDict = dict()
        if calcdValuesDict is None:
            calcdValuesDict = dict()
        outJSONFilePath = os.path.join(outputPath, outputFileName)
        jsonData = self.getJSONDictDefaultMetaData(
            productsStr, validMaskImage, footprintCalc, calcdValuesDict, outFilesDict
        )
        sensorInfo = jsonData["SensorInfo"]
        sensorInfo["Row"] = self.row
        sensorInfo["Path"] = self.path
        sensorInfo["SensorID"] = self.sensorID
        sensorInfo["SpacecraftID"] = self.spacecraftID
        acqDict = jsonData["AcquasitionInfo"]
        acqDict["EarthSunDistance"] = self.earthSunDistance
        imgInfo = dict()
        imgInfo["CloudCover"] = self.cloudCover
        imgInfo["CloudCoverLand"] = self.cloudCoverLand
        imgInfo["CellSizeRefl"] = self.gridCellSizeRefl
        imgInfo["CellSizeTherm"] = self.gridCellSizeTherm
        jsonData["ImageInfo"] = imgInfo

        with open(outJSONFilePath, "w") as outfile:
            json.dump(
                jsonData,
                outfile,
                sort_keys=True,
                indent=4,
                separators=(",", ": "),
                ensure_ascii=False,
            )

    def expectedImageDataPresent(self):
        imageDataPresent = True

        if not os.path.exists(self.band1File):
            imageDataPresent = False
        if not os.path.exists(self.band2File):
            imageDataPresent = False
        if not os.path.exists(self.band3File):
            imageDataPresent = False
        if not os.path.exists(self.band4File):
            imageDataPresent = False
        if not os.path.exists(self.band5File):
            imageDataPresent = False
        if not os.path.exists(self.band6File):
            imageDataPresent = False
        if not os.path.exists(self.band7File):
            imageDataPresent = False

        return imageDataPresent

    def hasThermal(self):
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
            "Landsat 5 TM does not provide any image masks, do not use the MASK option."
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
        raise ARCSIException("Image data does not need resampling")

    def sharpenLowResRadImgBands(self, inputImg, outputImage, outFormat):
        raise ARCSIException("Image sharpening is not available for this sensor.")

    def generateValidImageDataMask(
        self, outputPath, outputMaskName, viewAngleImg, outFormat
    ):
        print("Create the valid data mask")
        tmpBaseName = os.path.splitext(outputMaskName)[0]
        tmpValidPxlMsk = os.path.join(outputPath, tmpBaseName + "vldpxlmsk.kea")
        outputImage = os.path.join(outputPath, outputMaskName)
        inImages = [
            self.band1File,
            self.band2File,
            self.band3File,
            self.band4File,
            self.band5File,
            self.band7File,
            self.band6File,
        ]
        rsgislib.imageutils.gen_valid_mask(
            input_imgs=inImages,
            output_img=tmpValidPxlMsk,
            gdalformat="KEA",
            no_data_val=0.0,
        )
        rsgislib.rastergis.pop_rat_img_stats(tmpValidPxlMsk, True, False, True)
        # Check there is valid data
        ratDS = gdal.Open(tmpValidPxlMsk, gdal.GA_ReadOnly)
        Histogram = rat.readColumn(ratDS, "Histogram")
        ratDS = None
        if Histogram.shape[0] < 2:
            raise ARCSIException("There is no valid data in this image.")
        if not os.path.exists(viewAngleImg):
            print("Calculate Image Angles.")
            imgInfo = rios.fileinfo.ImageInfo(tmpValidPxlMsk)
            corners = fmask.landsatangles.findImgCorners(tmpValidPxlMsk, imgInfo)
            nadirLine = fmask.landsatangles.findNadirLine(corners)
            extentSunAngles = fmask.landsatangles.sunAnglesForExtent(
                imgInfo, self.fmaskMTLInfo
            )
            satAzimuth = fmask.landsatangles.satAzLeftRight(nadirLine)
            fmask.landsatangles.makeAnglesImage(
                tmpValidPxlMsk,
                viewAngleImg,
                nadirLine,
                extentSunAngles,
                satAzimuth,
                imgInfo,
            )
            dataset = gdal.Open(viewAngleImg, gdal.GA_Update)
            if not dataset is None:
                dataset.GetRasterBand(1).SetDescription("SatelliteAzimuth")
                dataset.GetRasterBand(2).SetDescription("SatelliteZenith")
                dataset.GetRasterBand(3).SetDescription("SolorAzimuth")
                dataset.GetRasterBand(4).SetDescription("SolorZenith")
            dataset = None
        rsgislib.imagecalc.band_math(
            outputImage,
            "(VA<14)&&(VM==1)?1:0",
            outFormat,
            rsgislib.TYPE_8UINT,
            [
                rsgislib.imagecalc.BandDefn("VA", viewAngleImg, 2),
                rsgislib.imagecalc.BandDefn("VM", tmpValidPxlMsk, 1),
            ],
        )

        rsgislib.imageutils.delete_gdal_layer(tmpValidPxlMsk)
        return outputImage

    def convertImageToRadiance(
        self, outputPath, outputReflName, outputThermalName, outFormat
    ):
        print("Converting to Radiance")
        outputReflImage = os.path.join(outputPath, outputReflName)
        outputThermalImage = None
        bandDefnSeq = list()

        lsBand = collections.namedtuple(
            "LSBand",
            [
                "band_name",
                "input_img",
                "img_band",
                "l_min",
                "l_max",
                "q_cal_min",
                "q_cal_max",
            ],
        )
        bandDefnSeq.append(
            lsBand(
                band_name="Blue",
                input_img=self.band1File,
                img_band=1,
                l_min=self.b1MinRad,
                l_max=self.b1MaxRad,
                q_cal_min=self.b1CalMin,
                q_cal_max=self.b1CalMax,
            )
        )
        bandDefnSeq.append(
            lsBand(
                band_name="Green",
                input_img=self.band2File,
                img_band=1,
                l_min=self.b2MinRad,
                l_max=self.b2MaxRad,
                q_cal_min=self.b2CalMin,
                q_cal_max=self.b2CalMax,
            )
        )
        bandDefnSeq.append(
            lsBand(
                band_name="Red",
                input_img=self.band3File,
                img_band=1,
                l_min=self.b3MinRad,
                l_max=self.b3MaxRad,
                q_cal_min=self.b3CalMin,
                q_cal_max=self.b3CalMax,
            )
        )
        bandDefnSeq.append(
            lsBand(
                band_name="NIR",
                input_img=self.band4File,
                img_band=1,
                l_min=self.b4MinRad,
                l_max=self.b4MaxRad,
                q_cal_min=self.b4CalMin,
                q_cal_max=self.b4CalMax,
            )
        )
        bandDefnSeq.append(
            lsBand(
                band_name="SWIR1",
                input_img=self.band5File,
                img_band=1,
                l_min=self.b5MinRad,
                l_max=self.b5MaxRad,
                q_cal_min=self.b5CalMin,
                q_cal_max=self.b5CalMax,
            )
        )
        bandDefnSeq.append(
            lsBand(
                band_name="SWIR2",
                input_img=self.band7File,
                img_band=1,
                l_min=self.b7MinRad,
                l_max=self.b7MaxRad,
                q_cal_min=self.b7CalMin,
                q_cal_max=self.b7CalMax,
            )
        )
        rsgislib.imagecalibration.landsat_to_radiance(
            outputReflImage, outFormat, bandDefnSeq
        )

        if not outputThermalName == None:
            outputThermalImage = os.path.join(outputPath, outputThermalName)
            bandDefnSeq = list()
            lsBand = collections.namedtuple(
                "LSBand",
                [
                    "band_name",
                    "input_img",
                    "img_band",
                    "l_min",
                    "l_max",
                    "q_cal_min",
                    "q_cal_max",
                ],
            )
            bandDefnSeq.append(
                lsBand(
                    band_name="ThermalB6",
                    input_img=self.band6File,
                    img_band=1,
                    l_min=self.b6MinRad,
                    l_max=self.b6MaxRad,
                    q_cal_min=self.b6CalMin,
                    q_cal_max=self.b6CalMax,
                )
            )
            rsgislib.imagecalibration.landsat_to_radiance(
                outputThermalImage, outFormat, bandDefnSeq
            )

        return outputReflImage, outputThermalImage

    def generateImageSaturationMask(self, outputPath, outputName, outFormat):
        print("Generate Saturation Image")
        outputImage = os.path.join(outputPath, outputName)

        lsBand = collections.namedtuple(
            "LSBand", ["band_name", "input_img", "img_band", "sat_val"]
        )
        bandDefnSeq = list()
        bandDefnSeq.append(
            lsBand(
                band_name="Blue",
                input_img=self.band1File,
                img_band=1,
                sat_val=self.b1CalMax,
            )
        )
        bandDefnSeq.append(
            lsBand(
                band_name="Green",
                input_img=self.band2File,
                img_band=1,
                sat_val=self.b2CalMax,
            )
        )
        bandDefnSeq.append(
            lsBand(
                band_name="Red",
                input_img=self.band3File,
                img_band=1,
                sat_val=self.b3CalMax,
            )
        )
        bandDefnSeq.append(
            lsBand(
                band_name="NIR",
                input_img=self.band4File,
                img_band=1,
                sat_val=self.b4CalMax,
            )
        )
        bandDefnSeq.append(
            lsBand(
                band_name="SWIR1",
                input_img=self.band5File,
                img_band=1,
                sat_val=self.b5CalMax,
            )
        )
        bandDefnSeq.append(
            lsBand(
                band_name="SWIR2",
                input_img=self.band7File,
                img_band=1,
                sat_val=self.b7CalMax,
            )
        )
        bandDefnSeq.append(
            lsBand(
                band_name="ThermalB6",
                input_img=self.band6File,
                img_band=1,
                sat_val=self.b6CalMax,
            )
        )

        rsgislib.imagecalibration.saturated_pixels_mask(
            outputImage, outFormat, bandDefnSeq
        )

        return outputImage

    def convertThermalToBrightness(
        self, inputRadImage, outputPath, outputName, outFormat, scaleFactor
    ):
        print("Converting to Thermal Brightness")
        outputThermalImage = os.path.join(outputPath, outputName)
        bandDefnSeq = list()

        lsBand = collections.namedtuple("LSBand", ["band_name", "img_band", "k1", "k2"])
        bandDefnSeq.append(
            lsBand(band_name="ThermalB6", img_band=1, k1=607.76, k2=1260.56)
        )
        rsgislib.imagecalibration.landsat_thermal_rad_to_brightness(
            inputRadImage,
            outputThermalImage,
            outFormat,
            rsgislib.TYPE_32INT,
            scaleFactor,
            bandDefnSeq,
        )
        return outputThermalImage

    def convertImageToTOARefl(
        self, inputRadImage, outputPath, outputName, outFormat, scaleFactor
    ):
        print("Converting to TOA")
        outputImage = os.path.join(outputPath, outputName)
        solarIrradianceVals = list()
        IrrVal = collections.namedtuple("SolarIrradiance", ["irradiance"])
        solarIrradianceVals.append(IrrVal(irradiance=1957.0))
        solarIrradianceVals.append(IrrVal(irradiance=1826.0))
        solarIrradianceVals.append(IrrVal(irradiance=1554.0))
        solarIrradianceVals.append(IrrVal(irradiance=1036.0))
        solarIrradianceVals.append(IrrVal(irradiance=215.0))
        solarIrradianceVals.append(IrrVal(irradiance=80.67))
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
            solarIrradianceVals,
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
        import rsgislib.imageutils

        try:
            outputImage = os.path.join(outputPath, outputCloudName)
            tmpBaseName = os.path.splitext(outputCloudName)[0]
            tmpBaseDIR = os.path.join(tmpPath, tmpBaseName)

            tmpDIRExisted = True
            if not os.path.exists(tmpBaseDIR):
                os.makedirs(tmpBaseDIR)
                tmpDIRExisted = False

            if (cloud_msk_methods is None) or (cloud_msk_methods == "FMASK"):
                tmpFMaskOut = os.path.join(tmpBaseDIR, tmpBaseName + "_pyfmaskout.kea")

                tmpThermalLayer = self.band6File
                if not rsgislib.imageutils.do_gdal_layers_have_same_proj(
                    inputThermalImage, self.band6File
                ):
                    tmpThermalLayer = os.path.join(
                        tmpBaseDIR, tmpBaseName + "_thermalresample.kea"
                    )
                    rsgislib.imageutils.resample_img_to_match(
                        inputThermalImage,
                        self.band6File,
                        tmpThermalLayer,
                        "KEA",
                        rsgislib.INTERP_CUBIC,
                        rsgislib.TYPE_32FLOAT,
                    )

                minCloudSize = 0
                cloudBufferDistance = 150
                shadowBufferDistance = 300

                fmaskFilenames = fmask.config.FmaskFilenames()
                fmaskFilenames.setTOAReflectanceFile(inputReflImage)
                fmaskFilenames.setThermalFile(tmpThermalLayer)
                fmaskFilenames.setSaturationMask(inputSatImage)
                fmaskFilenames.setOutputCloudMaskFile(tmpFMaskOut)

                thermalGain1040um = (self.b6MaxRad - self.b6MinRad) / (
                    self.b6CalMax - self.b6CalMin
                )
                thermalOffset1040um = self.b6MinRad - self.b6CalMin * thermalGain1040um
                thermalBand1040um = 0
                thermalInfo = fmask.config.ThermalFileInfo(
                    thermalBand1040um,
                    thermalGain1040um,
                    thermalOffset1040um,
                    607.76,
                    1260.56,
                )

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

                fmaskConfig = fmask.config.FmaskConfig(fmask.config.FMASK_LANDSAT47)
                fmaskConfig.setTOARefScaling(float(scaleFactor))
                fmaskConfig.setThermalInfo(thermalInfo)
                fmaskConfig.setAnglesInfo(anglesInfo)
                fmaskConfig.setKeepIntermediates(False)
                fmaskConfig.setVerbose(True)
                fmaskConfig.setTempDir(tmpBaseDIR)
                fmaskConfig.setMinCloudSize(minCloudSize)
                fmaskConfig.setEqn17CloudProbThresh(
                    fmask.config.FmaskConfig.Eqn17CloudProbThresh
                )
                fmaskConfig.setEqn20NirSnowThresh(
                    fmask.config.FmaskConfig.Eqn20NirSnowThresh
                )
                fmaskConfig.setEqn20GreenSnowThresh(
                    fmask.config.FmaskConfig.Eqn20GreenSnowThresh
                )

                # Work out a suitable buffer size, in pixels, dependent on the resolution of the input TOA image
                toaImgInfo = rios.fileinfo.ImageInfo(inputReflImage)
                fmaskConfig.setCloudBufferSize(
                    int(cloudBufferDistance / toaImgInfo.xRes)
                )
                fmaskConfig.setShadowBufferSize(
                    int(shadowBufferDistance / toaImgInfo.xRes)
                )

                fmask.fmask.doFmask(fmaskFilenames, fmaskConfig)

                rsgislib.imagecalc.image_math(
                    tmpFMaskOut,
                    outputImage,
                    "(b1==2)?1:(b1==3)?2:0",
                    outFormat,
                    rsgislib.TYPE_8UINT,
                )
            elif cloud_msk_methods == "LSMSK":
                if (self.bandQAFile == "") or (not os.path.exists(self.bandQAFile)):
                    raise ARCSIException(
                        "The QA band is not present - cannot use this for cloud masking."
                    )

                bqa_img_file = self.bandQAFile
                if not rsgislib.imageutils.do_gdal_layers_have_same_proj(
                    bqa_img_file, inputReflImage
                ):
                    bqa_img_file = os.path.join(tmpBaseDIR, tmpBaseName + "_BQA.kea")
                    rsgislib.imageutils.resample_img_to_match(
                        inputReflImage,
                        self.bandQAFile,
                        bqa_img_file,
                        "KEA",
                        rsgislib.INTERP_NEAREST_NEIGHBOUR,
                        rsgislib.TYPE_16UINT,
                        no_data_val=0,
                        multicore=False,
                    )

                if self.collection_num == 1:
                    exp = (
                        "(b1==752)||(b1==756)||(b1==760)||(b1==764)?1:"
                        "(b1==928)||(b1==932)||(b1==936)||(b1==940)||(b1==960)||(b1==964)||(b1==968)||(b1==972)?2:0"
                    )
                    rsgislib.imagecalc.image_math(
                        bqa_img_file, outputImage, exp, outFormat, rsgislib.TYPE_8UINT
                    )
                elif self.collection_num == 2:
                    import rsgislib.imagecalibration.sensorlvl2data

                    c2_bqa_ind_img_file = os.path.join(
                        tmpBaseDIR, tmpBaseName + "c2_qa_ind_bands.kea"
                    )
                    rsgislib.imagecalibration.sensorlvl2data.parse_landsat_c2_qa_pixel_img(
                        bqa_img_file, c2_bqa_ind_img_file, gdalformat="KEA"
                    )
                    band_defns = list()
                    band_defns.append(
                        rsgislib.imagecalc.BandDefn(
                            "DilatedCloud", c2_bqa_ind_img_file, 2
                        )
                    )
                    band_defns.append(
                        rsgislib.imagecalc.BandDefn("Cloud", c2_bqa_ind_img_file, 4)
                    )
                    band_defns.append(
                        rsgislib.imagecalc.BandDefn(
                            "CloudShadow", c2_bqa_ind_img_file, 5
                        )
                    )
                    rsgislib.imagecalc.band_math(
                        outputImage,
                        "(DilatedCloud == 1)||(Cloud == 1)?1:(CloudShadow == 1)?2:0",
                        "KEA",
                        rsgislib.TYPE_8UINT,
                        band_defns,
                    )
                else:
                    raise ARCSIException(
                        "Can only read Collection 1 and 2 cloud masks."
                    )

            else:
                raise ARCSIException(
                    "Landsat only has FMASK and LSMSK cloud masking options; option provided is unknown."
                )

            if outFormat == "KEA":
                rsgislib.rastergis.pop_rat_img_stats(outputImage, True, True)
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

            if not self.debugMode:
                if not tmpDIRExisted:
                    shutil.rmtree(tmpBaseDIR, ignore_errors=True)

            return outputImage, None
        except Exception as e:
            raise e

    def createCloudMaskDataArray(self, inImgDataArr):
        return inImgDataArr

    def defineDarkShadowImageBand(self):
        return 4

    def calc6SCoefficients(
        self, aeroProfile, atmosProfile, grdRefl, surfaceAltitude, aotVal, useBRDF
    ):
        sixsCoeffs = numpy.zeros((6, 6), dtype=numpy.float32)
        # Set up 6S model
        s = Py6S.SixS()
        s.atmos_profile = atmosProfile
        s.aero_profile = aeroProfile
        s.ground_reflectance = grdRefl
        s.geometry = Py6S.Geometry.Landsat_TM()
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

        # Band 1
        s.wavelength = Py6S.Wavelength(Py6S.PredefinedWavelengths.LANDSAT_TM_B1)
        s.run()
        sixsCoeffs[0, 0] = float(s.outputs.values["coef_xa"])
        sixsCoeffs[0, 1] = float(s.outputs.values["coef_xb"])
        sixsCoeffs[0, 2] = float(s.outputs.values["coef_xc"])
        sixsCoeffs[0, 3] = float(s.outputs.values["direct_solar_irradiance"])
        sixsCoeffs[0, 4] = float(s.outputs.values["diffuse_solar_irradiance"])
        sixsCoeffs[0, 5] = float(s.outputs.values["environmental_irradiance"])

        # Band 2
        s.wavelength = Py6S.Wavelength(Py6S.PredefinedWavelengths.LANDSAT_TM_B2)
        s.run()
        sixsCoeffs[1, 0] = float(s.outputs.values["coef_xa"])
        sixsCoeffs[1, 1] = float(s.outputs.values["coef_xb"])
        sixsCoeffs[1, 2] = float(s.outputs.values["coef_xc"])
        sixsCoeffs[1, 3] = float(s.outputs.values["direct_solar_irradiance"])
        sixsCoeffs[1, 4] = float(s.outputs.values["diffuse_solar_irradiance"])
        sixsCoeffs[1, 5] = float(s.outputs.values["environmental_irradiance"])

        # Band 3
        s.wavelength = Py6S.Wavelength(Py6S.PredefinedWavelengths.LANDSAT_TM_B3)
        s.run()
        sixsCoeffs[2, 0] = float(s.outputs.values["coef_xa"])
        sixsCoeffs[2, 1] = float(s.outputs.values["coef_xb"])
        sixsCoeffs[2, 2] = float(s.outputs.values["coef_xc"])
        sixsCoeffs[2, 3] = float(s.outputs.values["direct_solar_irradiance"])
        sixsCoeffs[2, 4] = float(s.outputs.values["diffuse_solar_irradiance"])
        sixsCoeffs[2, 5] = float(s.outputs.values["environmental_irradiance"])

        # Band 4
        s.wavelength = Py6S.Wavelength(Py6S.PredefinedWavelengths.LANDSAT_TM_B4)
        s.run()
        sixsCoeffs[3, 0] = float(s.outputs.values["coef_xa"])
        sixsCoeffs[3, 1] = float(s.outputs.values["coef_xb"])
        sixsCoeffs[3, 2] = float(s.outputs.values["coef_xc"])
        sixsCoeffs[3, 3] = float(s.outputs.values["direct_solar_irradiance"])
        sixsCoeffs[3, 4] = float(s.outputs.values["diffuse_solar_irradiance"])
        sixsCoeffs[3, 5] = float(s.outputs.values["environmental_irradiance"])

        # Band 5
        s.wavelength = Py6S.Wavelength(Py6S.PredefinedWavelengths.LANDSAT_TM_B5)
        s.run()
        sixsCoeffs[4, 0] = float(s.outputs.values["coef_xa"])
        sixsCoeffs[4, 1] = float(s.outputs.values["coef_xb"])
        sixsCoeffs[4, 2] = float(s.outputs.values["coef_xc"])
        sixsCoeffs[4, 3] = float(s.outputs.values["direct_solar_irradiance"])
        sixsCoeffs[4, 4] = float(s.outputs.values["diffuse_solar_irradiance"])
        sixsCoeffs[4, 5] = float(s.outputs.values["environmental_irradiance"])

        # Band 6
        s.wavelength = Py6S.Wavelength(Py6S.PredefinedWavelengths.LANDSAT_TM_B7)
        s.run()
        sixsCoeffs[5, 0] = float(s.outputs.values["coef_xa"])
        sixsCoeffs[5, 1] = float(s.outputs.values["coef_xb"])
        sixsCoeffs[5, 2] = float(s.outputs.values["coef_xc"])
        sixsCoeffs[5, 3] = float(s.outputs.values["direct_solar_irradiance"])
        sixsCoeffs[5, 4] = float(s.outputs.values["diffuse_solar_irradiance"])
        sixsCoeffs[5, 5] = float(s.outputs.values["environmental_irradiance"])

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

        return outputImage, elevAOTCoeffs

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
        s.geometry = Py6S.Geometry.Landsat_TM()
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
        s.wavelength = Py6S.Wavelength(Py6S.PredefinedWavelengths.LANDSAT_TM_B1)
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
        dataset.GetRasterBand(1).SetDescription("Blue")
        dataset.GetRasterBand(2).SetDescription("Green")
        dataset.GetRasterBand(3).SetDescription("Red")
        dataset.GetRasterBand(4).SetDescription("NIR")
        dataset.GetRasterBand(5).SetDescription("SWIR1")
        dataset.GetRasterBand(6).SetDescription("SWIR2")
        dataset = None

    def cleanLocalFollowProcessing(self):
        print("")
