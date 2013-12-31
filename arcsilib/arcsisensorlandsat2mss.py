"""
Module that contains the ARCSILandsat2MSSSensor class.
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
#           the pre-processing operations within ARCSI to the landsat 2 MSS
#           datasets.
#
# Author: Pete Bunting
# Email: pfb@aber.ac.uk
# Date: 13/07/2013
# Version: 1.0
#
# History:
# Version 1.0 - Created.
#
############################################################################

# import abstract base class stuff
from .arcsisensor import ARCSIAbstractSensor
# Import the ARCSI exception class
from .arcsiexception import ARCSIException
# Import the ARCSI utilities class
from .arcsiutils import ARCSIUtils
# Import the datetime module
import datetime
# Import the GDAL/OGR spatial reference library
from osgeo import osr
from osgeo import ogr
# Import OS path module for manipulating the file system 
import os.path
# Import the RSGISLib Image Calibration Module.
import rsgislib.imagecalibration
# Import the collections module
import collections
# Import the py6s module for running 6S from python.
import Py6S
# Import the python maths library
import math

class ARCSILandsat2MSSSensor (ARCSIAbstractSensor):
    """
    A class which represents the landsat 2 MSS sensor to read
    header parameters and apply data processing operations.
    """
    def __init__(self):
        ARCSIAbstractSensor.__init__(self)
        self.sensor = "LS2MSS"
        self.band4File = ""
        self.band5File = ""
        self.band6File = ""
        self.band7File = ""
        self.row = 0
        self.path = 0
        
        self.b4CalMin = 0
        self.b4CalMax = 0
        self.b5CalMin = 0
        self.b5CalMax = 0
        self.b6CalMin = 0
        self.b6CalMax = 0
        self.b7CalMin = 0
        self.b7CalMax = 0
        
        self.b4MinRad = 0.0
        self.b4MaxRad = 0.0
        self.b5MinRad = 0.0
        self.b5MaxRad = 0.0
        self.b6MinRad = 0.0
        self.b6MaxRad = 0.0
        self.b7MinRad = 0.0
        self.b7MaxRad = 0.0
    
    def extractHeaderParameters(self, inputHeader, wktStr):
        """
        Understands and parses the Landsat MTL header files
        """
        try:
            print("Reading header file")
            hFile = open(inputHeader, 'r')
            headerParams = dict()
            for line in hFile:
                line = line.strip()
                if line:
                    lineVals = line.split('=')
                    if len(lineVals) == 2:
                        if (lineVals[0].strip() != "GROUP") or (lineVals[0].strip() != "END_GROUP"):
                            headerParams[lineVals[0].strip()] = lineVals[1].strip().replace('"','')
            hFile.close()
            print("Extracting Header Values")
            # Get the sensor info.
            if (headerParams["SPACECRAFT_ID"] == "LANDSAT_2") and (headerParams["SENSOR_ID"] == "MSS"):
                self.sensor = "LS2MSS"
            else:
                raise ARCSIException("Do no recognise the spacecraft and sensor or combination.")
            
            # Get row/path
            self.row = int(headerParams["WRS_ROW"])
            self.path = int(headerParams["WRS_PATH"])
            
            # Get date and time of the acquisition
            acData = headerParams["DATE_ACQUIRED"].split('-')
            acTime = headerParams["SCENE_CENTER_TIME"].split(':')
            secsTime = acTime[2].split('.')
            self.acquisitionTime = datetime.datetime(int(acData[0]), int(acData[1]), int(acData[2]), int(acTime[0]), int(acTime[1]), int(secsTime[0]))
            
            self.solarZenith = 90-float(headerParams["SUN_ELEVATION"])
            self.solarAzimuth = float(headerParams["SUN_AZIMUTH"])
            
            # Get the geographic lat/long corners of the image.
            self.latTL = float(headerParams["CORNER_UL_LAT_PRODUCT"])
            self.lonTL = float(headerParams["CORNER_UL_LON_PRODUCT"])
            self.latTR = float(headerParams["CORNER_UR_LAT_PRODUCT"])
            self.lonTR = float(headerParams["CORNER_UR_LON_PRODUCT"])
            self.latBL = float(headerParams["CORNER_LL_LAT_PRODUCT"])
            self.lonBL = float(headerParams["CORNER_LL_LON_PRODUCT"])
            self.latBR = float(headerParams["CORNER_LR_LAT_PRODUCT"])
            self.lonBR = float(headerParams["CORNER_LR_LON_PRODUCT"])
            
            # Get the projected X/Y corners of the image
            self.xTL = float(headerParams["CORNER_UL_PROJECTION_X_PRODUCT"])
            self.yTL = float(headerParams["CORNER_UL_PROJECTION_Y_PRODUCT"])
            self.xTR = float(headerParams["CORNER_UR_PROJECTION_X_PRODUCT"])
            self.yTR = float(headerParams["CORNER_UR_PROJECTION_Y_PRODUCT"])
            self.xBL = float(headerParams["CORNER_LL_PROJECTION_X_PRODUCT"])
            self.yBL = float(headerParams["CORNER_LL_PROJECTION_Y_PRODUCT"])
            self.xBR = float(headerParams["CORNER_LR_PROJECTION_X_PRODUCT"])
            self.yBR = float(headerParams["CORNER_LR_PROJECTION_Y_PRODUCT"])
            
            # Get projection
            inProj = osr.SpatialReference()
            if (headerParams["MAP_PROJECTION"] == "UTM") and (headerParams["DATUM"] == "WGS84") and (headerParams["ELLIPSOID"] == "WGS84"):
                utmZone = int(headerParams["UTM_ZONE"])
                utmCode = "WGS84UTM" + str(utmZone) + str("N")
                #print("UTM: ", utmCode)
                inProj.ImportFromEPSG(self.epsgCodes[utmCode])
            else:
                raise ARCSIException("Expecting Landsat to be projected in UTM with datum=WGS84 and ellipsoid=WGS84.")
            
            # Check image is square!
            if not ((self.xTL == self.xBL) and (self.yTL == self.yTR) and (self.xTR == self.xBR) and (self.yBL == self.yBR)):
                raise ARCSIException("Image is not square in projected coordinates.")
            
            self.xCentre = self.xTL + ((self.xTR - self.xTL)/2)
            self.yCentre = self.yBR + ((self.yTL - self.yBR)/2)
            
            wgs84latlonProj = osr.SpatialReference()
            wgs84latlonProj.ImportFromEPSG(4326)
            
            wktPt = 'POINT(%s %s)' % (self.xCentre, self.yCentre)
            #print(wktPt)
            point = ogr.CreateGeometryFromWkt(wktPt)
            point.AssignSpatialReference(inProj)
            point.TransformTo(wgs84latlonProj)
            #print(point)
            
            self.latCentre = point.GetY()
            self.lonCentre = point.GetX()
            
            #print("Lat: " + str(self.latCentre) + " Long: " + str(self.lonCentre))
            
            filesDIR = os.path.dirname(inputHeader)
            
            self.band4File = os.path.join(filesDIR, headerParams["FILE_NAME_BAND_4"])
            self.band5File = os.path.join(filesDIR, headerParams["FILE_NAME_BAND_5"])
            self.band6File = os.path.join(filesDIR, headerParams["FILE_NAME_BAND_6"])
            self.band7File = os.path.join(filesDIR, headerParams["FILE_NAME_BAND_7"])
            
            self.b4CalMin = float(headerParams["QUANTIZE_CAL_MIN_BAND_4"])
            self.b4CalMax = float(headerParams["QUANTIZE_CAL_MAX_BAND_4"])
            self.b5CalMin = float(headerParams["QUANTIZE_CAL_MIN_BAND_5"])
            self.b5CalMax = float(headerParams["QUANTIZE_CAL_MAX_BAND_5"])
            self.b6CalMin = float(headerParams["QUANTIZE_CAL_MIN_BAND_6"])
            self.b6CalMax = float(headerParams["QUANTIZE_CAL_MAX_BAND_6"])
            self.b7CalMin = float(headerParams["QUANTIZE_CAL_MIN_BAND_7"])
            self.b7CalMax = float(headerParams["QUANTIZE_CAL_MAX_BAND_7"])
            
            self.b4MinRad = float(headerParams["RADIANCE_MINIMUM_BAND_4"])
            self.b4MaxRad = float(headerParams["RADIANCE_MAXIMUM_BAND_4"])
            self.b5MinRad = float(headerParams["RADIANCE_MINIMUM_BAND_5"])
            self.b5MaxRad = float(headerParams["RADIANCE_MAXIMUM_BAND_5"])
            self.b6MinRad = float(headerParams["RADIANCE_MINIMUM_BAND_6"])
            self.b6MaxRad = float(headerParams["RADIANCE_MAXIMUM_BAND_6"])
            self.b7MinRad = float(headerParams["RADIANCE_MINIMUM_BAND_7"])
            self.b7MaxRad = float(headerParams["RADIANCE_MAXIMUM_BAND_7"])
            
        except Exception as e:
            raise e
        
    def generateOutputBaseName(self):
        """
        Provides an implementation for the landsat sensor
        """
        rowpath = "r" + str(self.row) + "p" + str(self.path)
        outname = self.defaultGenBaseOutFileName()
        outname = outname + str("_") + rowpath
        return outname
        
    def convertImageToRadiance(self, outputPath, outputName, outFormat):
        print("Converting to Radiance")
        outputImage = os.path.join(outputPath, outputName)
        bandDefnSeq = list()
        lsBand = collections.namedtuple('LSBand', ['bandName', 'fileName', 'bandIndex', 'lMin', 'lMax', 'qCalMin', 'qCalMax'])
        bandDefnSeq.append(lsBand(bandName="Green", fileName=self.band4File, bandIndex=1, lMin=self.b4MinRad, lMax=self.b4MaxRad, qCalMin=self.b4CalMin, qCalMax=self.b4CalMax))
        bandDefnSeq.append(lsBand(bandName="Red", fileName=self.band5File, bandIndex=1, lMin=self.b5MinRad, lMax=self.b5MaxRad, qCalMin=self.b5CalMin, qCalMax=self.b5CalMax))
        bandDefnSeq.append(lsBand(bandName="NIR1", fileName=self.band6File, bandIndex=1, lMin=self.b6MinRad, lMax=self.b6MaxRad, qCalMin=self.b6CalMin, qCalMax=self.b6CalMax))
        bandDefnSeq.append(lsBand(bandName="NIR2", fileName=self.band7File, bandIndex=1, lMin=self.b7MinRad, lMax=self.b7MaxRad, qCalMin=self.b7CalMin, qCalMax=self.b7CalMax))
        rsgislib.imagecalibration.landsat2Radiance(outputImage, outFormat, bandDefnSeq)
        return outputImage
    
    def convertImageToTOARefl(self, inputRadImage, outputPath, outputName, outFormat):
        print("Converting to TOA")
        outputImage = os.path.join(outputPath, outputName)
        solarIrradianceVals = list()
        IrrVal = collections.namedtuple('SolarIrradiance', ['irradiance'])
        solarIrradianceVals.append(IrrVal(irradiance=1829.0))
        solarIrradianceVals.append(IrrVal(irradiance=1539.0))
        solarIrradianceVals.append(IrrVal(irradiance=1268.0))
        solarIrradianceVals.append(IrrVal(irradiance=886.6))
        rsgislib.imagecalibration.radiance2TOARefl(inputRadImage, outputImage, outFormat, rsgislib.TYPE_16UINT, 1000, self.acquisitionTime.year, self.acquisitionTime.month, self.acquisitionTime.day, self.solarZenith, solarIrradianceVals)
        return outputImage
        
    def convertImageToSurfaceReflSglParam(self, inputRadImage, outputPath, outputName, outFormat, aeroProfile, atmosProfile, grdRefl, surfaceAltitude, aotVal, useBRDF):
        print("Converting to Surface Reflectance")
        outputImage = os.path.join(outputPath, outputName)
        
        Band6S = collections.namedtuple('Band6SCoeff', ['band', 'aX', 'bX', 'cX'])
        imgBandCoeffs = list()
        # Set up 6S model
        s = Py6S.SixS()
        s.atmos_profile = atmosProfile
        s.aero_profile = aeroProfile
        s.ground_reflectance = grdRefl
        s.geometry = Py6S.Geometry.Landsat_TM()
        s.geometry.month = self.acquisitionTime.month
        s.geometry.day = self.acquisitionTime.day
        s.geometry.gmt_decimal_hour = float(self.acquisitionTime.hour) + float(self.acquisitionTime.minute)/60.0
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
        s.wavelength = Py6S.Wavelength(Py6S.SixSHelpers.PredefinedWavelengths.LANDSAT_MSS_B1)
        s.run()
        imgBandCoeffs.append(Band6S(band=1, aX=s.outputs.values['coef_xa'], bX=s.outputs.values['coef_xb'], cX=s.outputs.values['coef_xc']))
        
        # Band 2
        s.wavelength = Py6S.Wavelength(Py6S.SixSHelpers.PredefinedWavelengths.LANDSAT_MSS_B2)
        s.run()
        imgBandCoeffs.append(Band6S(band=2, aX=s.outputs.values['coef_xa'], bX=s.outputs.values['coef_xb'], cX=s.outputs.values['coef_xc']))
        
        # Band 3
        s.wavelength = Py6S.Wavelength(Py6S.SixSHelpers.PredefinedWavelengths.LANDSAT_MSS_B3)
        s.run()
        imgBandCoeffs.append(Band6S(band=3, aX=s.outputs.values['coef_xa'], bX=s.outputs.values['coef_xb'], cX=s.outputs.values['coef_xc']))
        
        # Band 4
        s.wavelength = Py6S.Wavelength(Py6S.SixSHelpers.PredefinedWavelengths.LANDSAT_MSS_B4)
        s.run()
        imgBandCoeffs.append(Band6S(band=4, aX=s.outputs.values['coef_xa'], bX=s.outputs.values['coef_xb'], cX=s.outputs.values['coef_xc']))
        
        for band in imgBandCoeffs:
            print(band)
        
        rsgislib.imagecalibration.apply6SCoeffSingleParam(inputRadImage, outputImage, outFormat, rsgislib.TYPE_16UINT, 1000, 0, True, imgBandCoeffs)
        
        return outputImage


    def estimateImageToAOD(self, inputTOAImage, outputPath, outputName, outFormat, tmpPath, aeroProfile, atmosProfile, grdRefl, surfaceAltitude, aotValMin, aotValMax):
        print("Not implemented\n")
        sys.exit()

