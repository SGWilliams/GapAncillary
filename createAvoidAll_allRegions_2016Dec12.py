##    createAvoidAll_allRegions.py
##    Curtis Belyea, Biologist / GIS Specialist; Biodiversity and Spatial Information Center;
##    USGS North Carolina Cooperative Fish and Wildlife Research Unit
##    Department of Applied Ecology, North Carolina State University
##
##
##    This script was written (as the second of two) to produce an avoidance mask used in modeling
##    predicted habitat for target vertabrate species.  The first produces regional road density raster
##    datasets used in this script.
##
##    In order to complete the avoid masks, this script does the following:
##
##        a.  Identifies groups of 5 cells or more in National Land Cover Dataset values 22, 23, 24 and expands (buffers)
##            them 120 meters (using the Euclidean Distance tool).
##        b.  Identifies groups of 5 cells or more in National Land Cover Dataset values 23, 24 and expands (buffers)
##            them 90 meters (using the Euclidean Distance tool).  
##        c.  Creates three new raster datasets that represent areas where the road density for the current
##            region is:
##                1.  greater than or equal to 25
##                2.  greater than or equal to 65
##                3.  greater than or equal to 130
##        d.  sets any cell location where a is true or c1 is true = 1.  Other cells = 0
##        e.  sets any cell location where b is true or c2 is true = 1.  Other cells = 0
##        f.  sets any cell location where b is true or c3 is true = 1.  Other cells = 0
##        g.  Sums up results of d, e and f to arrive at the final avoid mask for the region, with values ranging
##            from 0-3.
##
##        This process is repeated for each of the 6 regions in the continental US.
##
##
##    Script inputs:
##        Snapraster
##        National Land Cover Dataset
##        region boundary raster (from previous script)
##        region road density raster (from previous script)
##
##    Outputs:
##        Temporary:
##            tmp1_urb01 - result of con statement: if NLCD = 23, 24, 25 make it 1; otherwise 0
##            tmp1_urb02 - result of region group on tmp1_urb01; 8 cell neighborhood
##            tmp1_urb03 - result of zonal statistics (sum) on tmp1_urb02
##            tmp1_urb04 - result of con statement: if tmp1_urb03 > 5 make it 1; otherwise 0
##            tmp1_urb05 - result of euclidean distance on tmp1_urb04; max 120 meters
##            tmp1_urb06 - result of con(isnull()) on tmp1_urb05.  if > 0, make 1 otherwise 0
##
##            tmp2_urb01 - result of con statement: if NLCD = 23, 24 make it 1; otherwise 0
##            tmp2_urb02 - result of region group on tmp1_urb01; 8 cell neighborhood
##            tmp2_urb03 - result of zonal statistics (sum) on tmp1_urb02
##            tmp2_urb04 - result of con statement: if tmp1_urb03 > 5 make it 1; otherwise 0
##            tmp2_urb05 - result of euclidean distance on tmp1_urb04; max 90 meters
##            tmp2_urb06 - result of con(isnull()) on tmp1_urb05.  if > 0, make 1 otherwise 0
##
##            tmp1_rdd01 - result of con statement: if road density >= 25: make it 1; otherwise 0
##            tmp2_rdd01 - result of con statement: if road density >= 65: make it 1; otherwise 0
##            tmp2_rdd01 - result of con statement: if road density >= 130: make it 1; otherwise 0
##
##            tmp1_avd01 - result of sum of tmp1_urb06 and tmp1_rdd01 (representing cells in 120 meter buffer of groups of 5 or more cells
##                         with NLCD values 22, 23, 24 - OR - any cell with road density >= 25
##            tmp1_avd02 - result of con statement: if tmp1_avd01 > 0 make it 1; otherwise 0
##
##            tmp2_avd01 - result of sum of tmp2_urb06 and tmp2_rdd01 (representing cells in 90 meter buffer of groups of 5 or more cells
##                         with NLCD values 23, 24 - OR - any cell with road density >= 65
##            tmp2_avd02 - result of con statement: if tmp2_avd01 > 0 make it 1; otherwise 0
##
##            tmp3_avd01 - result of sum of tmp2_urb06 and tmp3_rdd01 (representing cells in 90 meter buffer of groups of 5 or more cells
##                         with NLCD values 23, 24 - OR - any cell with road density >= 130
##            tmp3_avd02 - result of con statement: if tmp3_avd01 > 0 make it 1; otherwise 0
##
##
##        Final:
##            i.e. "avoid_mask_ne" - final regional avoid mask: sum of tmp1_avd02 + tmp2_avd02 + tmp3_avd02
##
##
##    File directory structure:
##        region directory
##            intDens(dir) - rounded, integer line density rasters (from previous script)
##            mosaic(dir) - final road density raster (from previous script)
##            scratch(dir) - output location for "tmp" datasets
##            roads.gdb (geodatabase of roads for HUC 6 polygons that intersect regional boundary)
##            feature classes of manually removed roads features
##            regional boundary raster
##            feature class of regional HUC 6 polygons
##            output location for final avoid mask for region
##
####################################################################################################################################

import os, sys, arcgisscripting, arcpy
from arcpy.sa import *
##  Check out spatial analyst extension
arcpy.CheckOutExtension("Spatial")
##  allow outputs to be overwritten
arcpy.env.overwriteOutput = True

gp = arcgisscripting.create(9.3)
gp2 = arcgisscripting.create()

##regions = ["gp","ne","nw","se","sw","um"]


for region in regions:

    try:
        scratch = "Y:\\TIGER\\reg_" + region + "\\scratch"
        arcpy.env.extent = "Y:\\TIGER\\reg_" + region + "\\" + region + "_bnd"
        arcpy.env.mask = "Y:\\TIGER\\reg_" + region + "\\" + region + "_bnd"

        ##  set snap raster
        arcpy.env.snapRaster = "P:\\Proj3\\USGap\\Ancillary\\snapnlcd_2016"

        nlcd = "P:\\Proj7\\Landcover\\NLCD\\NLCD2011\\landcover\\nlcd_2011_landcover_2011__snapped_to_GAP.img"

        density = "Y:\\TIGER\\reg_" + region + "\\mosaic\\" + region + "_rddens"

        ##  set workspace
        arcpy.env.workspace = "Y:\\TIGER\\reg_" + region
        ## for low level avoid:

        tmp1_urb01 = Con((Raster(nlcd) > 21) & (Raster(nlcd) < 25),1,0)
        tmp1_urb01.save(scratch + "\\tmp1_urb01")

        tmp1_urb02 = RegionGroup(scratch + "\\tmp1_urb01","EIGHT","WITHIN","ADD_LINK","0")
        tmp1_urb02.save(scratch + "\\tmp1_urb02")

        tmp1_urb03 = ZonalStatistics(scratch + "\\tmp1_urb02","VALUE",scratch + "\\tmp1_urb01","SUM","DATA")
        tmp1_urb03.save(scratch + "\\tmp1_urb03")

        tmp1_urb04 = Con(Raster(scratch + "\\tmp1_urb03") > 5, 1)
        tmp1_urb04.save(scratch + "\\tmp1_urb04")

        tmp1_urb05 = EucDistance(scratch + "\\tmp1_urb04","120","30")
        tmp1_urb05.save(scratch + "\\tmp1_urb05")

        tmp1_urb06 = Con(IsNull(scratch + "\\tmp1_urb05"),0,1)
        tmp1_urb06.save(scratch + "\\tmp1_urb06")



        ## for med/high level avoid:

        tmp2_urb01 = Con((Raster(nlcd) > 22) & (Raster(nlcd) < 25),1,0)
        tmp2_urb01.save(scratch + "\\tmp2_urb01")

        tmp2_urb02 = RegionGroup(scratch + "\\tmp2_urb01","EIGHT","WITHIN","ADD_LINK","0")
        tmp2_urb02.save(scratch + "\\tmp2_urb02")

        tmp2_urb03 = ZonalStatistics(scratch + "\\tmp2_urb02","VALUE",scratch + "\\tmp2_urb01","SUM","DATA")
        tmp2_urb03.save(scratch + "\\tmp2_urb03")

        tmp2_urb04 = Con(Raster(scratch + "\\tmp2_urb03") > 5, 1)
        tmp2_urb04.save(scratch + "\\tmp2_urb04")

        tmp2_urb05 = EucDistance(scratch + "\\tmp2_urb04","90","30")
        tmp2_urb05.save(scratch + "\\tmp2_urb05")

        tmp2_urb06 = Con(IsNull(scratch + "\\tmp2_urb05"),0,1)
        tmp2_urb06.save(scratch + "\\tmp2_urb06")




        ##  reclass rd dens
        tmp1_rdd01 = Con(Raster(density) >= 25,1,0)
        tmp1_rdd01.save(scratch + "\\tmp1_rdd01")

        tmp2_rdd01 = Con(Raster(density) >= 65,1,0)
        tmp2_rdd01.save(scratch + "\\tmp2_rdd01")

        tmp3_rdd01 = Con(Raster(density) >= 130,1,0)
        tmp3_rdd01.save(scratch + "\\tmp3_rdd01")


        ##  add result of con(isnull(eucDist90) and low rd dens reclass
        ##  reclass to 1,0
        tmp1_avd01 = Plus(scratch + "\\tmp1_urb06",scratch + "\\tmp1_rdd01")
        tmp1_avd01.save(scratch + "\\tmp1_avd01")

        tmp1_avd02 = Con(Raster(scratch + "\\tmp1_avd01") > 0,1,0)
        tmp1_avd02.save(scratch + "\\tmp1_avd02")

        ##  add result of con(isnull(eucDist120) and med rd dens reclass
        ##  reclass to 1,0
        tmp2_avd01 = Plus(scratch + "\\tmp2_urb06",scratch + "\\tmp2_rdd01")
        tmp2_avd01.save(scratch + "\\tmp2_avd01")


        tmp2_avd02 = Con(Raster(scratch + "\\tmp2_avd01") > 0,1,0)
        tmp2_avd02.save(scratch + "\\tmp2_avd02")




        ##  add result of con(isnull(eucDist120) and high rd dens reclass
        ##  reclass to 1,0

        tmp3_avd01 = Plus(scratch + "\\tmp2_urb06",scratch + "\\tmp3_rdd01")
        tmp3_avd01.save(scratch + "\\tmp3_avd01")


        tmp3_avd02 = Con(Raster(scratch + "\\tmp3_avd01") > 0,1,0)
        tmp3_avd02.save(scratch + "\\tmp3_avd02")


        avdMask = Raster(scratch + "\\tmp1_avd02") + Raster(scratch + "\\tmp2_avd02") + Raster(scratch + "\\tmp3_avd02")
        avdMask.save("avoid_mask_" + region)


    except:
        print region + " failed."
