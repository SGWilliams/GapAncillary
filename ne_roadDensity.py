##    ne_roadDensity.py
##    Curtis Belyea, Biologist / GIS Specialist; Biodiversity and Spatial Information Center;
##    US FWS North Carolina Cooperative, North Carolina State University Department of Applied Ecology
##
##    This script was written as the first of two to produce an avoidance mask used in modeling
##    predicted habitat for target vertabrate species.  With this script, a regional road density
##    raster dataset is created for one region of the CONUS.  A separate script has been written for
##    each of the six regions of the CONUS due to the need for individualization for each region.  
##
##    Road density is calculated for each of the 6 digit Hydrologic Code (HUC 6 polygons)
##    which intersect the regional boundary.  Rounding and a Con(IsNull()) statement are then used to
##    produce a raster for each HUC6 which has zero values instead of NoData for the full extent and
##    area of the buffered region.  CellStatistics(Maximum) for all of those rasters then returns one
##    single raster for the region with integer road density values for all of the HUC 6 polygons.
##
##    Only roads with paved surfaces were to be considered in calculating road density for use in
##    producing the avoid mask, but no viable option exists for surface-type attribution.  As proxy
##    for surface type, due to the fact that TIGER roads attribution poorly adheres to their own
##    classification system at best, efforts were made to remove roads most likely to have dirt, gravel
##    or other non-paved surfaces through other attribution selections.  One important facet was excluding road
##    features with MTFCC code S1400 which also have no entry in the "FULLNAME" field.  Roads such as logging
##    roads, oil field roads and driveways often meet these two criteria.
##
##    In addition, attribute tables for the roads that intersect each huc were visually inspected for
##    names which indicate appropriateness for removal.  This would include roads which have MTFCC code
##    of S1400 and a value in "FULLNAME" such "Forest road".  Exclusions in such instances were employed
##    through attribution selection in the script on an individual basis for HUC 6 polygons and are commented.
##
##    Lastly, in a few places, roads included in the 2011 TIGER roads clearly did not exist in reality,
##    (the most recent Google Earth imagery, typically 2016 was used as verification for this).
##    In these cases, the offending road features were manually selected in ArcMap and deleted as no
##    viable option existed to exclude them through attribute selection.  Since their absence would easily
##    go unnoticed, prior to deletion these erroneous roads were copied to a new feature class dataset to provide
##    record of what was removed.  Feature classes were appropriately named (i.e. "130202_rdsManuallyRemoved.shp")
##    where the first six characters supply the corresponding 6 digit HUC code.
##
##    The second stage is executed through createAvoidAll_allRegions.py
##
##    Script inputs:
##        2011 US Census Bureau TIGER roads datasets which were appended regionally to all HUC polygons which
##        intersect the region polygon
##        HUC (Hydrologic Unit Code) 6 polygons
##        Snapraster
##        Region extent polygon
##
##    Outputs:
##        Temporary:
##            feature classes of roads intersecting each HUC 6 polygon
##            "raw" line density raster for each temp roads feature class: with an applied scale factor of 10,000.
##                The output units for these rasters are meters per 10,000 square meters.  With the methodology used,
##                the output units without a scale factor would be meters per square meters.  With the scale factor,
##                the output value becomes 215.42 meters per 10,000 square meter.  The scale factor was applied so that
##                upon conversion to an integer raster the result would have a more convenient range of values
##                used in creating density-threshold masks.
##            rounded, integer line density raster, with extent and area of region and NoData values set to 0
##                for each HUC
##            raster of region boundary
##        Final:
##            regional integer line density raster ("gp_rddens")
##
##
##    File directory structure:
##        region directory
##            intDens(dir) - rounded, integer line density rasters
##            mosaic(dir) - final line density raster (folder admittedly poorly named)
##            scratch(dir) - roads feature classes by huc as well as "raw" line density calculations
##            roads.gdb
##            feature classes of manually removed roads features
##            regional boundary raster
##            feature class of regional HUC 6 polygons
##
##
##  April 5, 2017 Updated script to utilize masks of buffered regions
##  before clipping results to actual region boundary
##
##  April 26, 2017 Updated documentation with explanation of line density analysis output units
##
#######################################################################################################

##  import modules
import sys, os, arcpy, datetime
from arcpy.sa import *
##  Check out spatial analyst extension
arcpy.CheckOutExtension("Spatial")
##  allow outputs to be overwritten
arcpy.env.overwriteOutput = True

##  define hucs
hucs = "Y:\\TIGER\\reg_ne\\ne_huc6rng_gap.shp"
##  define roads
rds = "Y:\\TIGER\\reg_ne\\roads.gdb\\ne_roads"
##  define scratch space
scratch = "Y:\\TIGER\\reg_ne\\scratch"
##  define results output directory
output = "Y:\\TIGER\\reg_ne\\intDens"

##  set snap raster
arcpy.env.snapRaster = "P:\\Proj3\\USGap\\Ancillary\\snapnlcd_2016"



##get list of hucs
hucIDs = []

##  create search cursor for huc layer
sc = arcpy.SearchCursor(hucs)
##  advance search cursor to first row
row = sc.next()

##  for each row
while row:
    ##  if the row's HUC 6 values is not in the list of hucIDs
    if row.getValue("HUC_6") not in hucIDs:
        ##  then add it
        hucIDs.append(row.getValue("HUC_6"))
    ##  advance the search cursor to the next row
    row = sc.next()

##  make feature layers of hucs and roads
arcpy.MakeFeatureLayer_management(hucs,"hucLayer")
arcpy.MakeFeatureLayer_management(rds,"rdsLayer")


print "Start time: " + str(datetime.datetime.now())
##  for each huc
for huc in hucIDs:
    try:
        print huc, str(hucIDs.index(huc)+ 1) + " / " + str(len(hucIDs))

        arcpy.env.mask = ""
        arcpy.env.extent = ""
        ##  select the huc with the current id
        arcpy.SelectLayerByAttribute_management("hucLayer","NEW_SELECTION","\"HUC_6\" = '" + huc + "'")

        ##  select roads by intersection with huc polygon and 575 meter distance
        arcpy.SelectLayerByLocation_management("rdsLayer","INTERSECT","hucLayer","575 Meters")

        ##  subset the selected roads to include only those categorized as paved, but excluding those with MTFCC code S1400 that do not have a name
        arcpy.SelectLayerByAttribute_management("rdsLayer", "SUBSET_SELECTION",'("MTFCC" = \'S1100\' OR "MTFCC" = \'S1200\' OR "MTFCC" = \'S1630\' OR "MTFCC" = \'S1730\') or ("MTFCC" = \'S1400\' AND "FULLNAME" <> \'\')')

        ##  remove problematic roads from a few hucs

        if huc == "020403":
            ##  remove some roads, easily identified from attributes as NUS Naval Surface Warfare Center roads
            arcpy.SelectLayerByAttribute_management("rdsLayer","REMOVE_FROM_SELECTION",'LOWER("FULLNAME") LIKE \'%dept of defense%\'')
            arcpy.MakeFeatureLayer_management("Y:\\TIGER\\reg_ne\\020403_rdsManuallyRemoved.shp","removeLayer")
            arcpy.SelectLayerByLocation_management("rdsLayer","ARE_IDENTICAL_TO","removeLayer","","REMOVE_FROM_SELECTION")
            arcpy.Delete_management("removeLayer")

        if huc == "050100":
            ##  remove some roads, easily identified from attributes as NUS Naval Surface Warfare Center roads
            arcpy.SelectLayerByAttribute_management("rdsLayer","REMOVE_FROM_SELECTION",'LOWER("FULLNAME")  LIKE \'%forest rd%\' OR LOWER("FULLNAME")  LIKE \'fr %\' OR LOWER("FULLNAME")  LIKE \'fs rd%\'')

        ##  copy selected roads to new feature class
        arcpy.CopyFeatures_management("rdsLayer",scratch + "\\" + huc + "tmp.shp")

        ##  dissolve new roads features to remove duplicates (there are more than you'd think)
        arcpy.Dissolve_management(scratch + "\\" + huc + "tmp.shp",scratch + "\\" + huc + "tmpDissolve.shp")

        ##  define line density expression
        expr = "linedensity(" + scratch + "\\" + huc + "tmpDissolve.shp" + ", none, 30, kernel, 10000, 564)"

        ##  execute line density through SOMA
        arcpy.gp.SingleOutputMapAlgebra_sa(expr, scratch + "\\Dens_" + str(huc))

    
    except:
        print huc, " failed."

##convert region polygon to raster
arcpy.FeatureToRaster_conversion("P:\\Proj3\\USGap\\Vert\\Model\\data\\regions\\regclip_ne.shp","FID","Y:\\TIGER\\reg_ne\\ne_bnd","30")

##  move workspace to directory with line density outputs
arcpy.env.workspace = scratch

##  set extent and mask to buffered region boundary grid
arcpy.env.mask = "Y:\\Avoid\\maskbuf_ne"
arcpy.env.extent = "Y:\\Avoid\\maskbuf_ne"

##  for each huc
for huc in hucIDs:
    ##  create a new grid with extent of region, adding 0.5 before converting to integer raster, to round appropriately
    ##  while simultaneously setting null values to 0    
    null = Con(IsNull(Int(Raster("dens_" + huc) + .5)),0,Int(Raster("dens_" + huc) + .5))
    ##  save new grid to output directory
    null.save(output + "\\null_" + huc)

##  move workspace to output directory
arcpy.env.workspace = output
##  create empty list of new rasters from previous step
nulls = []

##  for each raster in the workspace
for raster in arcpy.ListDatasets():
    ##  if it's name starts with "null_"
    if raster.startswith("null_"):
        ##  append it to the list
        nulls.append(raster)

##  use cell statistics to find max value of all null rasters at each pixel for region
arcpy.gp.cellstatistics_sa(nulls,"Y:\\TIGER\\reg_ne\\mosaic\\ne_rddens","MAXIMUM")

##  print stop time
print "Stop time: " + str(datetime.datetime.now())