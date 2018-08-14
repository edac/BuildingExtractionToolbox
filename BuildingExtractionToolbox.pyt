import arcpy
import os
import datetime
from arcpy import env
from time import sleep
from arcpy.sa import *
import glob

arcpy.env.overwriteOutput = True
timestamp = datetime.datetime.now()


class Toolbox(object):
    def __init__(self):
        self.label = "Building Toolbox"
        self.alias = "ArcGIS Building Toolbox"

        # List of tool classes associated with this toolbox
        self.tools = [Building_Extractor, NDVIBuilding_Filter, Building_Filter ]


class Building_Extractor(object):
    def __init__(self):
        self.label = "Step 1 - Building Object Extractor"
        self.description = "This tool creates raster building objects using LiDAR LAS 1.4 tiles.  The tool requires both the LAS and bare earth DEM tiles (as an Imagine .img, TIFF or ESRI Grid format) to be in separate folders and sharing the same naming convention.  This tool takes these files and a user defined minimum height above ground for a building to create the raster objects."
        self.canRunInBackground = False

    def getParameterInfo(self):

     # Input parameters
        lasdir = arcpy.Parameter(displayName="LAS Input Directory", name="lasdir",
                                 datatype="DEFolder", parameterType="Required", direction="Input")
        demdir = arcpy.Parameter(displayName="DEM Input Directory", name="demdir",
                                 datatype="DEFolder", parameterType="Required", direction="Input")
        outputdir = arcpy.Parameter(displayName="Output Directory", name="outputdir",
                                    datatype="DEFolder", parameterType="Required", direction="Input")
        spectral_detail = arcpy.Parameter(displayName="Spectral Detail", name="spectral_detail", datatype="GPDouble",
                                          parameterType="Required", direction="Input", category="Image Segmentation Parameters")
      # setting default value
        spectral_detail.value = 15.5
        spatial_detail = arcpy.Parameter(displayName="Spatial Detail", name="spatial_detail", datatype="GPLong",
                                         parameterType="Required", direction="Input", category="Image Segmentation Parameters")
      # setting default value
        spatial_detail.value = 15
        min_segment_size = arcpy.Parameter(displayName="Min Segment Size", name="min_segment_size", datatype="GPLong",
                                           parameterType="Required", direction="Input", category="Image Segmentation Parameters")
      # setting default value
        min_segment_size.value = 10
        height = arcpy.Parameter(displayName="Minimum Height (As Measured in Elevation Units)", name="height",
                                 datatype="GPDouble", parameterType="Required", direction="Input", category="Minimum Rooftop Height")
      # setting default value
        height.value = 2.0

        binningmethod = arcpy.Parameter(
            displayName="Sampling Method and Void Filling Method",
            name="binningmethod",
            datatype="GPString",
            parameterType="Required",
            direction="Input",
            category="LAS Dataset To Raster Parameters"
        )
        binningmethod.value = "BINNING MINIMUM NONE"
        binningmethod.filter.type = "ValueList"
        binningmethod.filter.list = ["BINNING AVERAGE NONE", "BINNING AVERAGE SIMPLE", "BINNING AVERAGE LINEAR", "BINNING AVERAGE NATURAL_NEIGHBOR", "BINNING MINIMUM NONE", "BINNING MINIMUM SIMPLE", "BINNING MINIMUM LINEAR", "BINNING MINIMUM NATURAL_NEIGHBOR", "BINNING MAXIMUM NONE",
                                     "BINNING MAXIMUM SIMPLE", "BINNING MAXIMUM LINEAR", "BINNING MAXIMUM NATURAL_NEIGHBOR", "BINNING IDW NONE", "BINNING IDW SIMPLE", "BINNING IDW LINEAR", "BINNING IDW NATURAL_NEIGHBOR", "BINNING NEAREST NONE", "BINNING NEAREST SIMPLE", "BINNING NEAREST LINEAR", "BINNING NEAREST NATURAL_NEIGHBOR"]

        lidarvalue = arcpy.Parameter(
            displayName="Data Attribute Used to Create DSM",
            name="lidarvalue",
            datatype="GPString",
            parameterType="Required",
            direction="Input",
            category="LAS Dataset To Raster Parameters"
        )
        lidarvalue.value = "ELEVATION"

        lidarvalue.filter.type = "ValueList"
        lidarvalue.filter.list = ["ELEVATION", "INTENSITY"]

        rasterouttype = arcpy.Parameter(
            displayName="DSM Raster Type",
            name="rasterouttype",
            datatype="GPString",
            parameterType="Required",
            direction="Input",
            category="LAS Dataset To Raster Parameters"
        )
        rasterouttype.value = "FLOAT"

        rasterouttype.filter.type = "ValueList"
        rasterouttype.filter.list = ["INT", "FLOAT"]

        samplingtype = arcpy.Parameter(
            displayName="Grid Cell Creation Method",
            name="samplingtype",
            datatype="GPString",
            parameterType="Required",
            direction="Input",
            category="LAS Dataset To Raster Parameters"
        )
        samplingtype.value = "CELLSIZE"

        samplingtype.filter.type = "ValueList"
        samplingtype.filter.list = ["OBSERVATIONS", "CELLSIZE"]

        samplingvalue = arcpy.Parameter(displayName="Grid Cell Spatial Resolution (As Measured in Horizontal Mapping Units)", name="samplingvalue", datatype="GPDouble",
                                        parameterType="Required", direction="Input", category="LAS Dataset To Raster Parameters")
        samplingvalue.value = 1

        parameters = [lasdir, demdir, outputdir, spectral_detail, spatial_detail, min_segment_size,
                      height, lidarvalue, binningmethod, rasterouttype, samplingtype, samplingvalue]
        return parameters

    def isLicensed(self):  # optional
        return True

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        arcpy.env.overwriteOutput = True
        arcpy.SetProgressor("default", "Working...", 0, 2, 1)
        env.workspace = arcpy.env.scratchFolder
        lasdir = parameters[0].valueAsText
        demdir = parameters[1].valueAsText
        outfolder = parameters[2].valueAsText
        spectral_detail = parameters[3].valueAsText
        spatial_detail = parameters[4].valueAsText
        min_segment_size = parameters[5].valueAsText
        band_indexes = ""
        height = parameters[6].valueAsText
        lidarval = parameters[7].valueAsText
        binningmethod = parameters[8].valueAsText
        data_type = parameters[9].valueAsText
        sampling_type = parameters[10].valueAsText
        sampling_value = parameters[11].valueAsText
        fulloutfolder = os.path.join(
            outfolder, timestamp.strftime('%Y%m%d%H%M%S'))

        arcpy.AddMessage("Creating output folder")
        if not os.path.exists(fulloutfolder):
            os.mkdir(fulloutfolder)

        files = [f for f in os.listdir(lasdir) if f.endswith(('.las', '.LAS'))]
        totalfiles = len(files)
        progress = 0
        for filename in files:  # os.listdir(lasdir):
            progress = progress+1
            arcpy.AddMessage("Running:" + filename + ", file " +
                             str(progress) + " of " + str(totalfiles))
            basename = filename.rstrip(".las")
            inputfile = os.path.join(lasdir, filename)
            outputfile = os.path.join(fulloutfolder, filename)
            arcpy.CreateLasDataset_management(
                inputfile, outputfile, create_las_prj="NO_FILES")
            lasLyr = arcpy.CreateUniqueName(basename)
            arcpy.management.MakeLasDatasetLayer(
                outputfile + "d", lasLyr, class_code=1, return_values='LAST RETURN')
            outimg = os.path.join(fulloutfolder, "lr"+basename+".img")
            arcpy.conversion.LasDatasetToRaster(
                lasLyr, outimg, lidarval, binningmethod, data_type, sampling_type, sampling_value, 1)
            arcpy.CheckOutExtension('Spatial')
            onlyfiles = glob.glob(os.path.join(demdir, basename+"*"))#[f for f in os.listdir(demdir) if os.path.isfile(os.path.join(demdir, f))]
            def firstvalid(filelist):
                for val in filelist:
                    if os.path.splitext(val)[1].lower()==".img":
                        return val
                    elif os.path.splitext(val)[1].lower()==".tif":
                        return val
                    elif os.path.splitext(val)[1].lower()==".tiff":
                        return val
                    elif os.path.splitext(val)[1].lower()==".grd":
                        return val
            validdem=firstvalid(onlyfiles)

            outMinus = Raster(outimg) - Raster(validdem)
            # outMinus.save(os.path.join(fulloutfolder,"lrdiff"+basename+".img"))
            outCon = Con(outMinus, outMinus, 0.00, "VALUE > " + height)
            # outCon.save(os.path.join(fulloutfolder,"lrdiffgt2"+basename+".img"))
            outSetNull = SetNull(outCon, outCon, "VALUE <= 0")
            outSetNull.save(os.path.join(
                fulloutfolder, "heightDSM"+basename+".img"))
            arcpy.AddMessage(
                "Segment Mean Shift phase. This will take some time. Be patient.")
            seg_raster = SegmentMeanShift(
                outSetNull, spectral_detail, spatial_detail,  min_segment_size)
            seg_raster.save(os.path.join(
                fulloutfolder, "lrdiffgt2is"+basename+".img"))
            CountoutCon = Con(seg_raster, seg_raster, 0.00, "COUNT < 10000")

            CountoutCon.save(os.path.join(fulloutfolder, "isobj"+basename+".img"))
            outPolygons = os.path.join(fulloutfolder, basename+".shp")
            field = "VALUE"
            arcpy.RasterToPolygon_conversion(
                CountoutCon, outPolygons, "NO_SIMPLIFY")
            outZonalStats = ZonalStatistics(
                outPolygons, "ID", outSetNull, "STD", "NODATA")
            outZonalStats.save(os.path.join(
                fulloutfolder, "isobjsd"+basename+".img"))
            arcpy.Delete_management(os.path.join(fulloutfolder, "isobj"+basename+".img"))
            arcpy.Delete_management(os.path.join(fulloutfolder, "lrdiffgt2is"+basename+".img"))
            arcpy.Delete_management(os.path.join(fulloutfolder, "lr"+basename+".img"))
            arcpy.Delete_management(os.path.join(fulloutfolder, basename+".shp"))
       
            arcpy.AddMessage("Finished:" + filename)
        return

        
class Building_Filter(object):
    def __init__(self):
        self.label = "Step 2A - SD Building Filter"
        self.description = "This tool creates building footprints from the output derived from the Building Object Extractor tool and then filters these objects using only the rooftop height standard deviation.  The results from these filters are then converted to vector polygons, filtered based on a minimum roof area and the edges cleaned up to create the building footprint polygons."
        self.canRunInBackground = False
    def getParameterInfo(self):

     # Input parameters
        inputdir = arcpy.Parameter(displayName="Input Directory", name="inputdir",
                                 datatype="DEFolder", parameterType="Required", direction="Input")
        outputdir = arcpy.Parameter(displayName="Output Directory", name="outputdir",
                                 datatype="DEFolder", parameterType="Required", direction="Input")
        threshold = arcpy.Parameter(displayName="Threshold", name="threshold",
                                 datatype="GPDouble", parameterType="Required", direction="Input")
        spectral_detail = arcpy.Parameter(displayName="Spectral Detail", name="spectral_detail", datatype="GPDouble",
                                          parameterType="Required", direction="Input", category="Image Segmentation Parameters")
      # setting default value
        spectral_detail.value = 15.5
        spatial_detail = arcpy.Parameter(displayName="Spatial Detail", name="spatial_detail", datatype="GPLong",
                                         parameterType="Required", direction="Input", category="Image Segmentation Parameters")
      # setting default value
        spatial_detail.value = 15
        min_segment_size = arcpy.Parameter(displayName="Min Segment Size", name="min_segment_size", datatype="GPLong",
                                           parameterType="Required", direction="Input", category="Image Segmentation Parameters")
        min_segment_size.value = 10
        # setting default value
        threshold.value = 1.5


        regularizationmethod = arcpy.Parameter(
            displayName="Regularization Method",
            name="regularizationmethod",
            datatype="GPString",
            parameterType="Required",
            direction="Input",
            category="Regularize Building Footprint"
        )
        regularizationmethod.value = "RIGHT_ANGLES"
        regularizationmethod.filter.type = "ValueList"
        regularizationmethod.filter.list = ["RIGHT_ANGLES","RIGHT_ANGLES_AND_DIAGONALS","ANY_ANGLE","CIRCLE"]

        tolerance = arcpy.Parameter(displayName="Tolerance", name="tolerance",
                                 datatype="GPDouble", parameterType="Required", direction="Input", category="Regularize Building Footprint")
        tolerance.value=2
        densification = arcpy.Parameter(displayName="Densification", name="densification",
                                 datatype="GPDouble", parameterType="Required", direction="Input", category="Regularize Building Footprint")
        densification.value=2

        shapearea = arcpy.Parameter(displayName="Minimum Building Footprint Area (in map units)", name="shapearea",
                                 datatype="GPLong", parameterType="Required", direction="Input", category="Footprint Size Threshold")
        shapearea.value=32


        parameters = [inputdir, outputdir, threshold, spectral_detail, spatial_detail, min_segment_size, regularizationmethod, tolerance, densification, shapearea]
        return parameters
    def isLicensed(self):  
        return True

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        arcpy.env.overwriteOutput = True
        arcpy.SetProgressor("default", "Working...", 0, 2, 1)
        env.workspace = arcpy.env.scratchFolder
        indir = parameters[0].valueAsText
        outfolder = parameters[1].valueAsText
        threshold = parameters[2].valueAsText
        spectral_detail = parameters[3].valueAsText
        spatial_detail = parameters[4].valueAsText
        min_segment_size = parameters[5].valueAsText
        regularizationmethod = parameters[6].valueAsText
        tolerance = parameters[7].valueAsText
        densification = parameters[8].valueAsText
        shapearea = parameters[9].valueAsText
        fulloutfolder = os.path.join(
            outfolder, timestamp.strftime('%Y%m%d%H%M%S'))

        arcpy.AddMessage("Creating output folder")
        if not os.path.exists(fulloutfolder):
           
                os.mkdir(fulloutfolder)
                
       

        files = [f for f in os.listdir(indir) if f.startswith('isobjsd') and f.endswith('img')]
        totalfiles = len(files)
        progress = 0
        for filename in files:  # os.listdir(lasdir):
            progress = progress+1
            arcpy.AddMessage("Running:" + filename + ", file " +
                             str(progress) + " of " + str(totalfiles))
            filepath= os.path.join(indir,filename)
            arcpy.AddMessage(filepath)
            basename = filename.rstrip(".img")
            arcpy.AddMessage(basename)
            outSetNull = SetNull(filepath, filepath, "VALUE >"+ threshold)
            outSetNull.save(os.path.join(fulloutfolder, "filt" + basename + ".img"))
            seg_raster = SegmentMeanShift((os.path.join(fulloutfolder, "filt" + basename + ".img")), spectral_detail, spatial_detail, min_segment_size)
            seg_raster.save((os.path.join(fulloutfolder, "seg" + basename + ".img")))
            CountoutCon = Con((os.path.join(fulloutfolder, "seg" + basename + ".img")), (os.path.join(fulloutfolder, "seg" + basename + ".img")), 0.00, "COUNT < 10000")
            outPolygons = os.path.join(fulloutfolder, basename+".shp")
            arcpy.RasterToPolygon_conversion(CountoutCon, outPolygons, "NO_SIMPLIFY")
            if not arcpy.Exists(os.path.join(fulloutfolder, "FinalBldgs.gdb")): 
                arcpy.CreateFileGDB_management(fulloutfolder, "FinalBldgs.gdb")
            finalbldgDissolve = os.path.join(fulloutfolder,"FinalBldgs.gdb","finalbldg"+basename+"Dissolve")
            RegBldgRightAngle = os.path.join(fulloutfolder,"FinalBldgs.gdb","RegBldg"+basename+"RightAngle")
            FinalBldg = os.path.join(fulloutfolder,"FinalBldgs.gdb","FinalBldg"+basename)
            arcpy.Dissolve_management(outPolygons, finalbldgDissolve, "", "", "SINGLE_PART", "DISSOLVE_LINES")
            arcpy.RegularizeBuildingFootprint_3d(finalbldgDissolve, RegBldgRightAngle, regularizationmethod, tolerance, densification, "0.25", "1.5", "0.1", "1000000")
            arcpy.Select_analysis(RegBldgRightAngle, FinalBldg, "Shape_Area >= "+ shapearea)



            #cleanup
            arcpy.Delete_management(os.path.join(fulloutfolder, "filt" + basename + ".img"))
            arcpy.Delete_management(os.path.join(fulloutfolder, "seg" + basename + ".img"))
            arcpy.Delete_management(os.path.join(fulloutfolder, basename+".shp"))
            arcpy.Delete_management(os.path.join(fulloutfolder,"FinalBldgs.gdb","finalbldg"+basename+"Dissolve"))
            arcpy.Delete_management(os.path.join(fulloutfolder,"FinalBldgs.gdb","RegBldg"+basename+"RightAngle"))
            




class NDVIBuilding_Filter(object):
    def __init__(self):
        self.label = "Step 2B - SD and NDVI Building Filter"
        self.description = "This tool creates building footprints from the output derived from the Building Object Extractor tool and then filters these objects using the rooftop height standard deviation and their average NDVI value.  The results from these filters are then converted to vector polygons, filtered based on a minimum roof area and the edges cleaned up to create the building footprint polygons."
        self.canRunInBackground = False
    def getParameterInfo(self):

     # Input parameters
        inputdir = arcpy.Parameter(displayName="Input Directory", name="inputdir",
                                 datatype="DEFolder", parameterType="Required", direction="Input")
        outputdir = arcpy.Parameter(displayName="Output Directory", name="outputdir",
                                 datatype="DEFolder", parameterType="Required", direction="Input")
        ndvifile = arcpy.Parameter(displayName="NDVI File", name="ndvifile",
                                 datatype="DEFile", parameterType="Required", direction="Input")

        threshold = arcpy.Parameter(displayName="Threshold", name="threshold",
                                 datatype="GPDouble", parameterType="Required", direction="Input")

        ndvithreshold = arcpy.Parameter(displayName=" NDVI Threshold", name="ndvithreshold",
                                 datatype="GPDouble", parameterType="Required", direction="Input")
        spectral_detail = arcpy.Parameter(displayName="Spectral Detail", name="spectral_detail", datatype="GPDouble",
                                          parameterType="Required", direction="Input", category="Image Segmentation Parameters")
      # setting default value
        spectral_detail.value = 15.5
        spatial_detail = arcpy.Parameter(displayName="Spatial Detail", name="spatial_detail", datatype="GPLong",
                                         parameterType="Required", direction="Input", category="Image Segmentation Parameters")
      # setting default value
        spatial_detail.value = 15
        min_segment_size = arcpy.Parameter(displayName="Min Segment Size", name="min_segment_size", datatype="GPLong",
                                           parameterType="Required", direction="Input", category="Image Segmentation Parameters")
        min_segment_size.value = 10
        # setting default value
        threshold.value = 1.5
        ndvithreshold.value = 105


        regularizationmethod = arcpy.Parameter(
            displayName="Regularization Method",
            name="regularizationmethod",
            datatype="GPString",
            parameterType="Required",
            direction="Input",
            category="Regularize Building Footprint"
        )
        regularizationmethod.value = "RIGHT_ANGLES"
        regularizationmethod.filter.type = "ValueList"
        regularizationmethod.filter.list = ["RIGHT_ANGLES","RIGHT_ANGLES_AND_DIAGONALS","ANY_ANGLE","CIRCLE"]

        tolerance = arcpy.Parameter(displayName="Tolerance", name="tolerance",
                                 datatype="GPDouble", parameterType="Required", direction="Input", category="Regularize Building Footprint")
        tolerance.value=2
        densification = arcpy.Parameter(displayName="Densification", name="densification",
                                 datatype="GPDouble", parameterType="Required", direction="Input", category="Regularize Building Footprint")
        densification.value=2

        shapearea = arcpy.Parameter(displayName="Minimum Building Footprint Area (in map units)", name="shapearea",
                                 datatype="GPLong", parameterType="Required", direction="Input", category="Footprint Size Threshold")
        shapearea.value=32


        parameters = [inputdir, outputdir, threshold, spectral_detail, spatial_detail, min_segment_size, ndvifile, ndvithreshold, regularizationmethod, tolerance, densification, shapearea ]
        return parameters
    def isLicensed(self):  # optional
        return True

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        arcpy.SetProgressor("default", "Working...", 0, 2, 1)
        env.workspace = arcpy.env.scratchFolder
        indir = parameters[0].valueAsText
        outfolder = parameters[1].valueAsText
        threshold = parameters[2].valueAsText
        
        spectral_detail = parameters[3].valueAsText
        spatial_detail = parameters[4].valueAsText
        min_segment_size = parameters[5].valueAsText
        ndvifile = parameters[6].valueAsText
        ndvithreshold = parameters[7].valueAsText
        regularizationmethod = parameters[8].valueAsText
        tolerance = parameters[9].valueAsText
        densification =  parameters[10].valueAsText
        shapearea =  parameters[11].valueAsText
        fulloutfolder = os.path.join(
            outfolder, timestamp.strftime('%Y%m%d%H%M%S'))

        arcpy.AddMessage("Creating output folder")
        if not os.path.exists(fulloutfolder):
           
                os.mkdir(fulloutfolder)

        files = [f for f in os.listdir(indir) if f.startswith('isobjsd') and f.endswith('img')]
        totalfiles = len(files)
        progress = 0
        for filename in files:  # os.listdir(lasdir):
            progress = progress+1
            arcpy.AddMessage("Running:" + filename + ", file " + str(progress) + " of " + str(totalfiles))
            filepath= os.path.join(indir,filename)
            arcpy.AddMessage(filepath)
            basename = filename.rstrip(".img")
            arcpy.AddMessage(basename)
            outSetNull = SetNull(filepath, filepath, "VALUE >"+ threshold)
            outSetNull.save(os.path.join(fulloutfolder, "filt" + basename + ".img"))
            seg_raster = SegmentMeanShift((os.path.join(fulloutfolder, "filt" + basename + ".img")), spectral_detail, spatial_detail, min_segment_size)
            seg_raster.save((os.path.join(fulloutfolder, "seg" + basename + ".img")))
            CountoutCon = Con((os.path.join(fulloutfolder, "seg" + basename + ".img")), (os.path.join(fulloutfolder, "seg" + basename + ".img")), 0.00, "COUNT < 10000")
            outPolygons = os.path.join(fulloutfolder, basename+".shp")
            arcpy.RasterToPolygon_conversion(CountoutCon, outPolygons, "NO_SIMPLIFY")


            outZonalStats = ZonalStatistics(outPolygons, "ID", ndvifile, "MEAN", "NODATA")
            outSetNull2 = SetNull(outZonalStats, outZonalStats, "VALUE > "+ ndvithreshold)
            outSetNull2.save(os.path.join(fulloutfolder, "NDVIfilt"+basename+".img"))
            seg_raster2 = SegmentMeanShift(os.path.join(fulloutfolder, "NDVIfilt"+basename+".img"), spectral_detail, spatial_detail, min_segment_size)
            seg_raster2.save((os.path.join(fulloutfolder, "seg" + basename + "2.img")))
            CountoutCon2 = Con(seg_raster2, seg_raster2, 0.00, "COUNT < 10000")
            outPolygons2 = os.path.join(fulloutfolder, basename+"-final.shp")
            arcpy.RasterToPolygon_conversion(CountoutCon2, outPolygons2, "NO_SIMPLIFY")
            if not arcpy.Exists(os.path.join(fulloutfolder, "FinalBldgs.gdb")):
                arcpy.CreateFileGDB_management(fulloutfolder, "FinalBldgs.gdb")
            finalbldgDissolve = os.path.join(fulloutfolder,"FinalBldgs.gdb","NDVIfinalbldg"+basename+"Dissolve")
            RegBldgRightAngle = os.path.join(fulloutfolder,"FinalBldgs.gdb","NDVIRegBldg"+basename+"RightAngle")
            FinalBldg = os.path.join(fulloutfolder,"FinalBldgs.gdb","NDVIFinalBldg"+basename)
            arcpy.Dissolve_management(outPolygons2, finalbldgDissolve, "", "", "SINGLE_PART", "DISSOLVE_LINES")
            arcpy.RegularizeBuildingFootprint_3d(finalbldgDissolve, RegBldgRightAngle, regularizationmethod, tolerance, densification, "0.25", "1.5", "0.1", "1000000")
            arcpy.Select_analysis(RegBldgRightAngle, FinalBldg, "Shape_Area >= " + shapearea)

            #cleanup
            arcpy.Delete_management(os.path.join(fulloutfolder, "filt" + basename + ".img"))
            arcpy.Delete_management(os.path.join(fulloutfolder, "seg" + basename + ".img"))
            arcpy.Delete_management(os.path.join(fulloutfolder, basename+".shp"))
            arcpy.Delete_management(os.path.join(fulloutfolder, "NDVIfilt"+basename+".img"))
            arcpy.Delete_management(os.path.join(fulloutfolder,"FinalBldgs.gdb","NDVIfinalbldg"+basename+"Dissolve"))
            arcpy.Delete_management(os.path.join(fulloutfolder,"FinalBldgs.gdb","NDVIRegBldg"+basename+"RightAngle"))
