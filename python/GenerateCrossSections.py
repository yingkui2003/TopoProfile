#-------------------------------------------------------------------------------
# Name: Paleo Ice Reconstruction with ice boundary
# Purpose: This tool calculates the ice thickness on flowlines based on the Excel flowline model introduced by the Benn and Houlton (2010)
# and the python model, GLaRe, developed by Pellitero et al.(2016). However, both of the above models require manually assign shear stress
# and F factors. It requires a lot of time and efforts to conduct a paleo ice reconstrcution work.
# 
# The new developed tool will automatically adjust shear stress and F factors based on the DEM and target features. For shear stress, the tool
# assumes one value for the whole glacier based on the recommendation of many previous studies. The shear stress is first derived based on a
# revised shearstress code from VOLTA based on the elevation distrbution of ice surface and then adjusted to reach the best fit to
# the ice boundary.
#
# This tool also automatically adjust the F factor (shape factor) based on the cross sections along the flowlines. There are two options to derive
# the F factor: one is based on the cross section area, ice-contact perimeter, and ice thickness; the other is based on the fit of the polynomial
# function introduced by Li et al (2012). This tool uses a maximum width to prevent the error long cross sections in some sections when
#  the direction is not paralell the valley direction and where the tributary valley joins the main valley. This tool also applied EU
# allocation method to make sure the cross section does not extend to the tributaries (May need to check if it is necessary because the width is
# already used to constrain the exent of the cross section).
# 
# Author:      Yingkui Li
# Created:     06/01/2021 to 06/02/2021
# Updated:     02/16/2023
# Copyright:   (c) Yingkui Li 2020
# Department of Geography, University of Tennessee
# Knoxville, TN 37996, USA
#-------------------------------------------------------------------------------

from __future__ import division
import arcpy
from arcpy import env
from arcpy.sa import *
import numpy as np
import math
import os, sys
import scipy
from scipy.spatial import cKDTree as KDTree
from scipy import ndimage
import arcpy.cartography as CA
import matplotlib.pyplot as plt

arcpy.env.overwriteOutput = True
arcpy.env.XYTolerance= "0.01 Meters"

##Check the python version to determine ArcGIS or ArcGIS Pro
ArcGISPro = 0
arcpy.AddMessage("The current python version is: " + str(sys.version_info[0]))
if sys.version_info[0] == 2:  ##For ArcGIS 10, need to check the 3D and Spatial Extensions
    try:
        if arcpy.CheckExtension("Spatial")=="Available":
            arcpy.CheckOutExtension("Spatial")
        else:
            raise Exception ("not extension available")
    except:
        raise Exception ("unable to check out extension")

    try:
        if arcpy.CheckExtension("3D")=="Available":
            arcpy.CheckOutExtension("3D")
        else:
            raise Exception ("not extension available")
    except:
        raise Exception ("unable to check out extension")
elif sys.version_info[0] == 3:  ##For ArcGIS Pro
    ArcGISPro = 1
else:
    raise Exception("Must be using Python 2.x or 3.x")
    exit()   

temp_workspace = "in_memory"  
if ArcGISPro:
    temp_workspace = "memory"



#------------------------------------------------------------------------------------------------------------
# This function check each line in the line feature and make sure the line is from low elevation to high
# elevation (Glacier streamline needs from low to high elevation in order to reconstruct the paleo ice thickness).
# It is revised from the codes by Pellitero et al.(2016) in GlaRe.
#------------------------------------------------------------------------------------------------------------
def Check_If_Flip_Line_Direction(line, dem):
    cellsize = arcpy.GetRasterProperties_management(dem,"CELLSIZEX")
    cellsize_int = int(float(cellsize.getOutput(0)))

    line3d = arcpy.env.scratchGDB + "\\line3d"
    arcpy.AddField_management(line, "Flip", "Long", "", "", "", "", "", "", "")

    arcpy.InterpolateShape_3d(dem, line, line3d, cellsize_int*3) 

    flip_list = []
    i = 0
    with arcpy.da.SearchCursor(line3d,["Shape@"]) as cursor:
        for row in cursor:
            startZ = row[0].firstPoint.Z
            endZ = row[0].lastPoint.Z

            if startZ >= endZ:  ##Flip = True use equal in case the start and end point are the same
                flip_list.append(1)
            else:  ##Flip = False
                flip_list.append(0)
            i += 1 

    del cursor
    if i>0:
        del row

    if sum(flip_list) > 0:
        with arcpy.da.UpdateCursor(line,["Flip"]) as cursor:
            i = 0
            for row in cursor:
                row[0] = flip_list[i]
                cursor.updateRow(row)
                i += 1 
        del row, cursor

        arcpy.MakeFeatureLayer_management(line, "lyrLines")
        arcpy.SelectLayerByAttribute_management("lyrLines", "NEW_SELECTION", '"Flip" > 0')
        #arcpy.AddMessage("The number of fliped lines is: " + str(sum(flip_list)))
        arcpy.FlipLine_edit("lyrLines")  ##Need to change to lyrLines
        arcpy.SelectLayerByAttribute_management("lyrLines", "CLEAR_SELECTION")

    arcpy.DeleteField_management (line, "Flip")
    arcpy.Delete_management(line3d)


#---------------------------------------------------------------------------------------------------------------
# This function calculates the distance between two points
#--------------------------------------------------------------------------------------------------------------- 
def Dist(x1,y1,x2,y2):
    return math.sqrt(math.pow(math.fabs(x1-x2),2)+math.pow(math.fabs(y1-y2),2))

# normalized K_curve Regression
def k_curve_fit(x, y):
    popt, pcov = curve_fit(k_curve, x, y)
    c = popt[0]
    # r-squared
    yhat = k_curve(x, c)             # or [p(z) for z in x]
    ybar = np.sum(y)/len(y)          # or sum(y)/len(y)
    ssreg = np.sum((yhat-ybar)**2)   # or sum([ (yihat - ybar)**2 for yihat in yhat])
    sstot = np.sum((y - ybar)**2)    # or sum([ (yi - ybar)**2 for yi in y])

    Residual = y - yhat
    R2 = ssreg / sstot
    if R2 > 1:
        R2 = 1/R2
    return (c, R2, Residual)

def turning_points_RDP(streamLength, streamZ): ###, turning_points = 10, cluster_radius = 200):
    ##use power law to fit the cross section and then calculate the different between Z and Z (curve)

    #LengthfromStart.reverse()
    max_Z = max(streamZ)
    min_Z = min(streamZ)
    max_len = streamLength[-1]
    min_len = streamLength[0]
    normalH = np.array([(y - min_Z)/(max_Z - min_Z) for y in streamZ]) ##+ 0.01
    normalLen = np.array([(y-min_len) /(max_len-min_len) for y in streamLength]) ##+ 0.01

    #plt.plot( normalLen, normalH)
    #plt.show()
    try: 
        fit_results = k_curve_fit(normalLen, normalH)
        c = fit_results[0]
        R2 = fit_results[1]
        residual = fit_results[2]
    except:
        #arcpy.AddMessage("Error K-curve")
        c = 100
    #arcpy.AddMessage(c)

    if c < 0:
        #plt.plot( streamLength, residual)
        #plt.show()
        max_idx = np.argmax(residual)
        #arcpy.AddMessage(max_idx)
        max_value = residual[max_idx] * (max_Z - min_Z) ##/linedist
        #arcpy.AddMessage(max_value)
    else:
        points = np.concatenate([np.expand_dims(streamLength, axis=-1), np.expand_dims(streamZ, axis=-1)], axis=-1)
        start = np.tile(np.expand_dims(points[0], axis=0), (points.shape[0], 1))
        end = np.tile(np.expand_dims(points[-1], axis=0), (points.shape[0], 1))
        #linedist = Dist(start[0][0],start[0][1],end[0][0],end[0][1])
        dist_point_to_line = np.cross(end - start, points - start, axis=-1) / np.linalg.norm(end - start, axis=-1)

        max_idx = np.argmax(dist_point_to_line)
        #arcpy.AddMessage(max_idx)
        max_value = dist_point_to_line[max_idx]##/linedist

    return max_idx, max_value

#------------------------------------------------------------------------------------------------------------
# This fuction creates perpendicular lines along the flowlines and then create a set of points along these
# perpendicular lines. This tool aslo extracts the elevation of each points. The max_width is used to limit
# the extent of these points. The division between the perp lines in different valleys are determined by the
# EU Allocation tool in ArcGIS.
#------------------------------------------------------------------------------------------------------------
def create_cross_sections(flowlinepoints, flowline, beddem, constrainboundary, eraseAreas, cellsize_float, half_width, spacing): ##, min_width, b_divide): 

    spatialref=arcpy.Describe(flowlinepoints).spatialReference
    width = half_width  ##only use the half the width for each side of the flowline

    ##Step 1: data preparation
    arcpy.AddMessage("Step 1: data preparation...")
    flowlinepointscp = arcpy.CopyFeatures_management(flowlinepoints, temp_workspace + "\\flowlinepointscp")
    arcpy.AddField_management(flowlinepointscp, 'PntID', 'Long', 6) ##OFID is not FID, but it is related to FID
    arcpy.CalculateField_management(flowlinepointscp,"PntID",str("!"+str(arcpy.Describe(flowlinepointscp).OIDFieldName)+"!"),"PYTHON_9.3")

    arcpy.AddField_management(flowline, 'SegmentID', 'Long', 6) 
    arcpy.CalculateField_management(flowline,"SegmentID",str("!"+str(arcpy.Describe(flowline).OIDFieldName)+"!"),"PYTHON_9.3")

    flowlinepointscopy = temp_workspace + "\\flowlinepointscopy"
    arcpy.SpatialJoin_analysis(flowlinepointscp, flowline, flowlinepointscopy, "JOIN_ONE_TO_ONE", "KEEP_COMMON", '#', "INTERSECT", "1 Meters", "#")

    arr=arcpy.da.FeatureClassToNumPyArray(flowlinepointscopy, ('SHAPE@X', 'SHAPE@Y', 'PntID', 'OFID', 'SegmentID'))
    segment_ids = np.array([item[4] for item in arr])
    unique_segment_ids = np.unique(segment_ids)

    ##Step 2: create perpendicular lines along the flowpoints
    arcpy.AddMessage("Step 2: create cross sections along the streamlines...")
    distance = width
        
    perpendiculars = arcpy.CreateFeatureclass_management(temp_workspace, "perpendiculars","POLYLINE", "","","", spatialref)
    arcpy.AddField_management(perpendiculars, 'FlowlineID', 'Long', 6) 
    arcpy.AddField_management(perpendiculars, 'SegmentID', 'Long', 6) 
    arcpy.AddField_management(perpendiculars, 'FlowPntID', 'Long', 6) 
    
    exclude_PntID = []
    new_line_cursor = arcpy.da.InsertCursor(perpendiculars, ('SHAPE@', 'FlowPntID' , 'FlowlineID', 'SegmentID'))
    
    for row in range(len(arr)):
        if row > 0: ##not the first point
            startx = arr[row][0]
            starty = arr[row][1]
            creat_perp = 1
            if row <(len(arr)-1):
                if arr[row][3]==arr[row+1][3] and arr[row][3]==arr[row-1][3]: ##inside points
                    endx = arr[row+1][0]
                    endy = arr[row+1][1]
                elif arr[row][3]==arr[row+1][3]: ##start tributary points
                    exclude_PntID.append(arr[row][2])
                    creat_perp = 0
                else: ##end of the tributary
                    endx = arr[row-1][0]
                    endy = arr[row-1][1]
            else:
                endx = arr[row-1][0]
                endy = arr[row-1][1]

            if creat_perp > 0:
                if starty==endy or startx==endx:
                    if starty == endy:
                        y1 = starty + distance
                        y2 = starty - distance
                        x1 = startx
                        x2 = startx
                    if startx == endx:
                        y1 = starty
                        y2 = starty 
                        x1 = startx + distance
                        x2 = startx - distance     
                else:
                    ##This is revised codes, much simplified!! 11/09/2024
                    offsety = distance * (startx - endx) / spacing
                    offsetx = distance * (starty - endy) / spacing
                    y1 = starty + offsety 
                    y2 = starty - offsety
                    x1 = startx - offsetx
                    x2 = startx + offsetx
 
                array = arcpy.Array([arcpy.Point(x1,y1),arcpy.Point(x2, y2)])
                polyline = arcpy.Polyline(array)
                lineID = arr[row][3]
                pntID = arr[row][2]
                segID = arr[row][4]
                new_line_cursor.insertRow([polyline, pntID, lineID, segID])
                
    del new_line_cursor
    '''
    ##Step 3: Use EU Alloation to cut perp lines
    arcpy.AddMessage("Step 3: Use EU Alloation to cut perp lines")
    ##divide watershed by different flowlines
    ##Delete the exclude_points
    if len(exclude_PntID) > 0:
        with arcpy.da.UpdateCursor(flowlinepointscopy, 'PntID') as cursor:
            for row in cursor:
                if row in exclude_PntID:
                    cursor.deleteRow()
        del row, cursor

    # Process: Euclidean Allocation
    eucAllocate = EucAllocation(flowlinepointscopy, "", "", cellsize_float,'SegmentID', "", "")
    
    arcpy.RasterToPolygon_conversion(eucAllocate, temp_workspace + "\\eucAllocate", "SIMPLIFY", "VALUE")

    ##Clip the cross section using the watershed boundary
    if constrainboundary != "":
        arcpy.Clip_analysis (temp_workspace + "\\eucAllocate", constrainboundary, temp_workspace + "\\clipAllocation")
        arcpy.CopyFeatures_management(temp_workspace + "\\clipAllocation", temp_workspace + "\\eucAllocate")

    selAllocation = temp_workspace + "\\selAllocation"
    sel_perp = temp_workspace + "\\sel_perp"
    clipedlines = temp_workspace + "\\clipedlines"

    for segment_id in unique_segment_ids:
        query = "gridcode = " + str(segment_id) 
        arcpy.Select_analysis (temp_workspace + "\\eucAllocate", selAllocation, query)
        query = "SegmentID = "+ str(segment_id) 
        arcpy.Select_analysis (perpendiculars, sel_perp, query)

        arcpy.Clip_analysis (sel_perp, selAllocation, temp_workspace + "\\clip_perp")

        if segment_id == unique_segment_ids[0]:
            arcpy.CopyFeatures_management(temp_workspace + "\\clip_perp", clipedlines)
        else:
            arcpy.Append_management(temp_workspace + "\\clip_perp", clipedlines, "NO_TEST" )
        

    
    ##Multipart to singleparts
    singlepartlines = arcpy.MultipartToSinglepart_management(clipedlines, temp_workspace + "\\singlepartlines")
    '''

    ##Step 3: Use EU Alloation to cut perp lines
    arcpy.AddMessage("Step 3: Cut and clean up cross sections...")
    ##divide watershed by different flowlines
    ##Delete the exclude_points
    if len(exclude_PntID) > 0:
        with arcpy.da.UpdateCursor(flowlinepointscopy, 'PntID') as cursor:
            for row in cursor:
                if row in exclude_PntID:
                    cursor.deleteRow()
        del row, cursor

    #create minimum bounding geometry, convex hull method"
    if constrainboundary != "":
        mbg = arcpy.MinimumBoundingGeometry_management(constrainboundary, temp_workspace + "\\mbg", "ENVELOPE", "ALL", "","NO_MBG_FIELDS")
    else:
        mbg = arcpy.MinimumBoundingGeometry_management(flowline, temp_workspace + "\\mbg", "ENVELOPE", "ALL", "","NO_MBG_FIELDS")
    ##create a 300 m buffer around the mbg
    buffer_dist = str(half_width) + " Meter"
    mbg_buf = arcpy.Buffer_analysis(mbg, temp_workspace + "\\mbg_buf", buffer_dist, "", "", "ALL")

    extDEM = ExtractByMask(beddem,mbg_buf)

    arcpy.env.extent = extDEM
    arcpy.env.cellSize = extDEM
    arcpy.env.snapRaster = extDEM ##setup snap raster

    fillDEM =Fill(extDEM)  ##Fill the sink first
    fdir = FlowDirection(fillDEM,"NORMAL") ##Flow direction
    facc = FlowAccumulation(fdir) ##Flow accmulation

    ##covert the flowlines to raster
    streamlink = temp_workspace + "\\streamlink"
    arcpy.conversion.FeatureToRaster(flowline, 'SegmentID', streamlink)
    outWs = Watershed(fdir, streamlink)
    arcpy.RasterToPolygon_conversion(outWs, temp_workspace + "\\eucAllocate", "SIMPLIFY", "VALUE")

    #arcpy.CopyFeatures_management(temp_workspace + "\\eucAllocate", "d:\\tempLyk\\eucAllocate.shp") 
    #arcpy.CopyFeatures_management(perpendiculars, "d:\\tempLyk\\perpendiculars.shp") 

    ##Clip the cross section using the watershed boundary
    if constrainboundary != "":
        arcpy.Clip_analysis (temp_workspace + "\\eucAllocate", constrainboundary, temp_workspace + "\\clipAllocation")
        arcpy.CopyFeatures_management(temp_workspace + "\\clipAllocation", temp_workspace + "\\eucAllocate")

    selAllocation = temp_workspace + "\\selAllocation"
    sel_perp = temp_workspace + "\\sel_perp"
    clipedlines = temp_workspace + "\\clipedlines"

    for segment_id in unique_segment_ids:
        query = "gridcode = " + str(segment_id) 
        arcpy.Select_analysis (temp_workspace + "\\eucAllocate", selAllocation, query)
        query = "SegmentID = "+ str(segment_id) 
        arcpy.Select_analysis (perpendiculars, sel_perp, query)

        arcpy.Clip_analysis (sel_perp, selAllocation, temp_workspace + "\\clip_perp")

        if segment_id == unique_segment_ids[0]:
            arcpy.CopyFeatures_management(temp_workspace + "\\clip_perp", clipedlines)
        else:
            arcpy.Append_management(temp_workspace + "\\clip_perp", clipedlines, "NO_TEST" )
        
    ##Multipart to singleparts
    singlepartline = arcpy.MultipartToSinglepart_management(clipedlines, temp_workspace + "\\singlepartline")
    fieldmappings = arcpy.FieldMappings()
    fieldmappings.addTable(singlepartline)
    singlepartlines = temp_workspace + "\\singlepartlines"
    arcpy.SpatialJoin_analysis(singlepartline, flowline, singlepartlines, "JOIN_ONE_TO_ONE", "KEEP_COMMON", fieldmappings, "INTERSECT", "1 Meters", "#")

    ##Step 4: Clean perp lines, so that only the part intersect the flowpoints are preserved
    #arcpy.AddMessage("Step 4: Clean perp lines, so that only the part intersect the flowpoints are preserved")
    lines = []
    flowlineids = []
    flowpntids = []
    keep = []
    with arcpy.da.SearchCursor(singlepartlines, ['SHAPE@', 'FlowlineID', 'FlowPntID']) as perps:
        i = 0
        for perp in perps:
            lines.append(perp[0])
            flowlineids.append(perp[1]) 
            flowpntids.append(perp[2])
            keep.append(0)
            i += 1
    if i > 0:
        del perp
    del perps

    flowpntid_arr = np.array(flowpntids)

    with arcpy.da.SearchCursor(flowlinepointscopy, ['SHAPE@', 'OFID', 'OID@']) as pnts:
        for pnt in pnts:
            index = np.where(flowpntid_arr == pnt[2])
            indexlen = len(index[0])
            if indexlen > 0:
                if indexlen > 1: ##more than two lines for a point, only select the lines intersected with the point
                    for a in range (indexlen):
                        value = index[0][a]
                        within = pnt[0].within(lines[value])
                        if within==True: 
                            keep[value] = 1
                            break
                else:
                    value = index[0][0]
                    keep[value] = 1
    del pnts, pnt

    with arcpy.da.UpdateCursor(singlepartlines, 'SHAPE@') as perps:
        i = 0
        for perp in perps:
            if keep[i] == 0:
                perps.deleteRow()
            i += 1
    if i>0:
        del perp
    del perps

    if eraseAreas != "":
        fieldmappings = arcpy.FieldMappings()
        fieldmappings.addTable(temp_workspace + "\\singlepartlines")
        arcpy.SpatialJoin_analysis(singlepartlines, eraseAreas, temp_workspace + "\\sel_perps", "JOIN_ONE_TO_ONE", "KEEP_COMMON", fieldmappings, "INTERSECT", "1 Meters", "#")
        arcpy.analysis.Erase(singlepartlines, temp_workspace + "\\sel_perps", temp_workspace + "\\final_perps")
        return temp_workspace + "\\final_perps"
    else:
        return singlepartlines

## This is the new function that nees to write
#------------------------------------------------------------------------------------------------------------
# This fuction is the whole process to reconstruct paleoice based on DEM, input flowlines, ice boundary, and default shear stress
#------------------------------------------------------------------------------------------------------------
def CreateCrossSections(BedDEM, inputflowline, constrainboundary, eraseAreas, spacing, half_width, AdjustProfile, min_width, min_height, b_divide, out_cross_sections, OutputConvexPoints):

    GlacierID = "GlacierID" ##This is an ID field in inputflowline to identify the flowline(s) for each glacier (maybe connected with multiple flowlines)

    icebndpolys = temp_workspace + "\\icebndpolys"

    #arcpy.CopyFeatures_management(iceboundary, icebndpolys)

    ####Flow direction and accumulation analysis
    cellsize = arcpy.GetRasterProperties_management(BedDEM,"CELLSIZEX")
    cellsize_float = float(cellsize.getOutput(0)) # use float cell size
    spatialref=arcpy.Describe(inputflowline).spatialReference #get spat ref from input    

    #Check and flip lines if necessary    
    arcpy.AddMessage("Checking flowline direction...")
    Check_If_Flip_Line_Direction(inputflowline, BedDEM)

    #Make a copy of the flowline
    flowlines = temp_workspace + "\\flowlines"
    flowline = temp_workspace + "\\flowline"
    flowline3dpoints = temp_workspace + "\\flowline3dpoints"
    selflowline3dpoints = temp_workspace + "\\selflowline3dpoints"
    singepoint = temp_workspace + "\\singepoint"
    ws = temp_workspace + "\\ws"
    wsflowpoints = temp_workspace + "\\wsflowpoints"
    icepolys = temp_workspace + "\\icepolys"
    icepolyselect = temp_workspace + "\\icepolyselect"
    mainflowline = temp_workspace + "\\mainflowline"
    flowlinecopy = temp_workspace + "\\flowlinecopy"
    allicepolys = temp_workspace + "\\allicepolys"
    singeflowpnts = temp_workspace + "\\singeflowpnts"
    icewatersheds = temp_workspace + "\\icewatersheds"

    arcpy.CopyFeatures_management(inputflowline, flowlines)

    arcpy.env.extent = BedDEM

    ###The process to ordering the flowlines
    #Obtain the height info for the start of each flowline
    height=[]
    lineid = []
    #glaciers = []
    with arcpy.da.SearchCursor(flowlines, ['SHAPE@']) as cursor:
        i = 0
        for row in cursor:
            Startpoint = row[0].firstPoint
            coord= str(Startpoint.X)+" "+str(Startpoint.Y)
            Cellvalue=arcpy.GetCellValue_management(BedDEM, coord)
            Startpoint.Z=Cellvalue.getOutput(0)
            height.append(Startpoint.Z)
            lineid.append(i)
            #glaciers.append(row[1])
            i += 1

    del cursor, row

    ##Order the line geometries in the list
    arcpy.AddMessage("Ordering flowlines...")
    arcpy.AddField_management(flowlines,"ProcessID","LONG",6)

    order = sorted(range(len(height)), key=lambda k: height[k])  ##order is the ID

    with arcpy.da.UpdateCursor(flowlines, "ProcessID") as cursor: ##Fix the assigning order issue
        i = 0
        for row in cursor:
            row[0] = order.index(i)
            cursor.updateRow(row)
            i += 1
    del row, cursor
    
    OriginalFID = []
    Points = []
    processPID = []
    #glaIds = []
    p = 0

    geometry = arcpy.CopyFeatures_management(flowlines, arcpy.Geometry()) 

    for i in order:
        Length = geometry[i].length
        Length = int(Length)
        try:
            rlist = xrange(0, Length, spacing)
        except: ##python 3 xrange is replaced by range
            rlist = range(0, Length, spacing)
        
        for j in rlist:            
            Points.append(geometry[i].positionAlongLine(j))
            OriginalFID.append(i)
            #glaIds.append(glaciers[i])
            processPID.append(p)
        p += 1

    PointsFile = temp_workspace + "\\PointsFile"
    arcpy.CopyFeatures_management(Points, PointsFile)

    ##Extract Values and xy coordinates from Elevation Raster
    ExtractValuesToPoints(PointsFile, BedDEM, flowline3dpoints,"INTERPOLATE", "VALUE_ONLY")

    ##Add original FID Field
    arcpy.AddField_management(flowline3dpoints, 'OFID', 'Long', 6) ##OFID is not FID, but it is related to FID
    arcpy.AddField_management(flowline3dpoints, 'ProcessID', 'Long', 6) ##OFID is not FID, but it is related to FID
    #arcpy.AddField_management(flowline3dpoints, GlacierID, 'Long', 6) ##OFID is not FID, but it is related to FID
    PointsCursor = arcpy.UpdateCursor(flowline3dpoints, ['OFID','ProcessID'])
    i = 0
    for row in PointsCursor:
        row.OFID = OriginalFID[i]
        row.ProcessID = processPID[i]
        #row.GlacierID = glaIds[i]
        PointsCursor.updateRow(row)
        i+=1
    del row, PointsCursor  
    
    arcpy.DeleteField_management(flowline3dpoints,["POINT_Z","POINT_M"])

    '''
    ##remove the sections without big regional difference within the thieson polygons
    ThiessenPolygons = arcpy.CreateThiessenPolygons_analysis(flowline3dpoints, temp_workspace + "\\ThiessenPolygons", "ALL")
    buffer_dist = "1500 Meter"
    big_buf = arcpy.Buffer_analysis(flowlines, temp_workspace + "\\big_buf", buffer_dist, "", "", "ALL")

    buffer_dist = "100 Meter"
    small_buf = arcpy.Buffer_analysis(flowlines, temp_workspace + "\\small_buf", buffer_dist, "", "", "ALL")

    ClippedBigThiessenPolygons = arcpy.Clip_analysis(ThiessenPolygons, big_buf, temp_workspace + "\\ClippedBigThiessenPolygons", "")
    arcpy.AddField_management(ClippedBigThiessenPolygons, 'PolyID', 'Long', 6)
    arcpy.CalculateField_management(ClippedBigThiessenPolygons,"PolyID",str("!"+str(arcpy.Describe(ClippedBigThiessenPolygons).OIDFieldName)+"!"),"PYTHON_9.3")
    outZSaT = ZonalStatisticsAsTable(ClippedBigThiessenPolygons, 'PolyID', Raster(BedDEM), temp_workspace + "\\zonalSATbig", "#", "ALL")
    #fieldList = ["PCT90"]
    #arcpy.JoinField_management(ClippedBigThiessenPolygons, 'PolyID', temp_workspace + "\\zonalSAT", 'PolyID', fieldList)

    ClippedSmallThiessenPolygons = arcpy.Clip_analysis(ThiessenPolygons, small_buf, temp_workspace + "\\ClippedSmallThiessenPolygons", "")
    arcpy.AddField_management(ClippedSmallThiessenPolygons, 'PolyID', 'Long', 6)
    arcpy.CalculateField_management(ClippedSmallThiessenPolygons,"PolyID",str("!"+str(arcpy.Describe(ClippedSmallThiessenPolygons).OIDFieldName)+"!"),"PYTHON_9.3")
    outZSaT = ZonalStatisticsAsTable(ClippedSmallThiessenPolygons, 'PolyID', Raster(BedDEM), temp_workspace + "\\zonalSATsmall", "#", "ALL")
    fieldList = ["PCT90"]
    arcpy.JoinField_management(ClippedSmallThiessenPolygons, 'PolyID', temp_workspace + "\\zonalSATsmall", 'PolyID', fieldList)
    arcpy.JoinField_management(ClippedSmallThiessenPolygons, 'PolyID', temp_workspace + "\\zonalSATbig", 'PolyID', fieldList)


    
    arcpy.CopyFeatures_management(ClippedSmallThiessenPolygons, "d:\\tempLyk\\ClippedSmallThiessenPolygons.shp")
    #arcpy.CopyFeatures_management(ClippedBigThiessenPolygons, "d:\\tempLyk\\ClippedBigThiessenPolygons.shp")
    '''
    
    ##Step 3: Use EU Alloation to cut perp lines
    #arcpy.AddMessage("Step 3: Cut perp lines using subwatersheds")
    ##divide watershed by different flowlines
    ##Delete the exclude_points
    '''
    if len(exclude_PntID) > 0:
        with arcpy.da.UpdateCursor(flowlinepointscopy, 'PntID') as cursor:
            for row in cursor:
                if row in exclude_PntID:
                    cursor.deleteRow()
        del row, cursor
    '''
    '''
    #create minimum bounding geometry, convex hull method"
    if constrainboundary != "":
        mbg = arcpy.MinimumBoundingGeometry_management(constrainboundary, temp_workspace + "\\mbg", "ENVELOPE", "ALL", "","NO_MBG_FIELDS")
    else:
        mbg = arcpy.MinimumBoundingGeometry_management(flowlines, temp_workspace + "\\mbg", "ENVELOPE", "ALL", "","NO_MBG_FIELDS")
    ##create a 300 m buffer around the mbg
    buffer_dist = str(width) + " Meter"
    mbg_buf = arcpy.Buffer_analysis(mbg, temp_workspace + "\\mbg_buf", buffer_dist, "", "", "ALL")

    extDEM = ExtractByMask(BedDEM,mbg_buf)

    arcpy.env.extent = extDEM
    arcpy.env.cellSize = extDEM
    arcpy.env.snapRaster = extDEM ##setup snap raster

    fillDEM =Fill(extDEM)  ##Fill the sink first
    fdir = FlowDirection(fillDEM,"NORMAL") ##Flow direction
    facc = FlowAccumulation(fdir) ##Flow accmulation

    ##covert the flowlines to raster
    streamlink = temp_workspace + "\\streamlink"
    arcpy.conversion.FeatureToRaster(flowlines, "ProcessID", streamlink)
    outWs = Watershed(fdir, streamlink)
    arcpy.RasterToPolygon_conversion(outWs, temp_workspace + "\\eucAllocate", "SIMPLIFY", "VALUE")

    arcpy.CopyFeatures_management(temp_workspace + "\\eucAllocate", "d:\\tempLyk\\eucAllocate.shp") 

    ##Clip the cross section using the watershed boundary
    if constrainboundary != "":
        arcpy.Clip_analysis (temp_workspace + "\\eucAllocate", constrainboundary, temp_workspace + "\\clipAllocation")
        arcpy.CopyFeatures_management(temp_workspace + "\\clipAllocation", temp_workspace + "\\eucAllocate")
    '''
    '''
    selAllocation = temp_workspace + "\\selAllocation"
    sel_perp = temp_workspace + "\\sel_perp"
    clipedlines = temp_workspace + "\\clipedlines"
    for segment_id in unique_segment_ids:
        query = "gridcode = " + str(segment_id) 
        arcpy.Select_analysis (temp_workspace + "\\eucAllocate", selAllocation, query)
        query = "SegmentID = "+ str(segment_id) 
        arcpy.Select_analysis (perpendiculars, sel_perp, query)

        arcpy.Clip_analysis (sel_perp, selAllocation, temp_workspace + "\\clip_perp")

        if segment_id == unique_segment_ids[0]:
            arcpy.CopyFeatures_management(temp_workspace + "\\clip_perp", clipedlines)
        else:
            arcpy.Append_management(temp_workspace + "\\clip_perp", clipedlines, "NO_TEST" )
    ##Multipart to singleparts
    singlepartline = arcpy.MultipartToSinglepart_management(clipedlines, temp_workspace + "\\singlepartline")
    fieldmappings = arcpy.FieldMappings()
    fieldmappings.addTable(singlepartline)
    singlepartlines = temp_workspace + "\\singlepartlines"
    arcpy.SpatialJoin_analysis(singlepartline, flowline, singlepartlines, "JOIN_ONE_TO_ONE", "KEEP_COMMON", fieldmappings, "INTERSECT", "1 Meters", "#")

    ##Step 4: Clean perp lines, so that only the part intersect the flowpoints are preserved
    arcpy.AddMessage("Step 4: Clean perp lines, so that only the part intersect the flowpoints are preserved")
    lines = []
    flowlineids = []
    flowpntids = []
    keep = []
    with arcpy.da.SearchCursor(singlepartlines, ['SHAPE@', 'FlowlineID', 'FlowPntID']) as perps:
        i = 0
        for perp in perps:
            lines.append(perp[0])
            flowlineids.append(perp[1]) 
            flowpntids.append(perp[2])
            keep.append(0)
            i += 1
    if i > 0:
        del perp
    del perps

    flowpntid_arr = np.array(flowpntids)

    with arcpy.da.SearchCursor(flowlinepointscopy, ['SHAPE@', 'OFID', 'OID@']) as pnts:
        for pnt in pnts:
            index = np.where(flowpntid_arr == pnt[2])
            indexlen = len(index[0])
            if indexlen > 0:
                if indexlen > 1: ##more than two lines for a point, only select the lines intersected with the point
                    for a in range (indexlen):
                        value = index[0][a]
                        within = pnt[0].within(lines[value])
                        if within==True: 
                            keep[value] = 1
                            break
                else:
                    value = index[0][0]
                    keep[value] = 1
    del pnts, pnt

    with arcpy.da.UpdateCursor(singlepartlines, 'SHAPE@') as perps:
        i = 0
        for perp in perps:
            if keep[i] == 0:
                perps.deleteRow()
            i += 1
    if i>0:
        del perp
    del perps

    '''


    exist_fields = [f.name for f in arcpy.ListFields(flowlines)] #List of current field names in outline layer
    line_fields = ["line_id"]
    for field in line_fields:
        if field not in exist_fields:
            arcpy.AddField_management(flowlines, field, "LONG")
    arcpy.CalculateField_management(flowlines,"line_id",str("!"+str(arcpy.Describe(flowlines).OIDFieldName)+"!"),"PYTHON_9.3")

    singlepartlines = create_cross_sections(flowline3dpoints, flowlines, BedDEM, constrainboundary, eraseAreas, cellsize_float, half_width, spacing)##, min_width, b_divide)

    ##refine the cross sections
    #arcpy.AddMessage("Step 5: Constrain Cross section widths...")
    final_cross_sections = temp_workspace + "\\final_cross_sections"
    if len(AdjustProfile) > 10:  ##Refine profiles 

        profile3D = temp_workspace + "\\profile3D"
        arcpy.InterpolateShape_3d(BedDEM, singlepartlines, profile3D)

        lowest_X_coord = []
        lowest_Y_coord = []
        with arcpy.da.SearchCursor(profile3D, ["SHAPE@"]) as cursor:
            for row in cursor: ##Loop for each line
                PointX = []
                PointY = []
                PointZ = []
                for part in row[0]:
                    for pnt in part:
                        if pnt:
                            PointX.append(pnt.X)
                            PointY.append(pnt.Y)
                            PointZ.append(pnt.Z)

                pointZArr = np.array(PointZ).astype(int)
                pntXarr = np.array(PointX)
                pntYarr = np.array(PointY)

                ##Get the X Y coordinates of the lowest point
                min_Z = min(pointZArr)
                lowest_X_coord.append (pntXarr[pointZArr == min_Z][0])
                lowest_Y_coord.append (pntYarr[pointZArr == min_Z][0])

        '''
        lowest_points = arcpy.CreateFeatureclass_management(temp_workspace, "lowest_points","POINT", "","","", spatialref)
        arcpy.AddField_management(lowest_points, 'PntID', 'Long', 6) 

        new_point_cursor = arcpy.da.InsertCursor(lowest_points, ('SHAPE@', 'PntID'))
        for i in range(len(lowest_X_coord)):
            pnt = arcpy.Point(lowest_X_coord[i],lowest_Y_coord[i])
            new_point_cursor.insertRow([pnt, i])
        del new_point_cursor        
        '''

        if "convex" in AdjustProfile:
            arcpy.AddMessage("Step 4: Cut the cross sections by the convex points on each side...")
        else:
            arcpy.AddMessage("Step 4: Cut the cross sections by the highest points on each side...")
        #FID_list = []
        X_coord = []
        Y_coord = []
        pntType = []
        FID = []
        Height = []
        side = []
        Length = []
        with arcpy.da.SearchCursor(temp_workspace + "\\profile3D", ["SHAPE@", "OID@"]) as cursor:
            for row in cursor: ##Loop for each line
                PointX = []
                PointY = []
                LengthfromStart = []
                PointZ = []
                oid = row[1]
                cumLength = 0
                for part in row[0]:
                    pntCount = 0
                    cumLength = 0
                    segmentdis = 0
                    for pnt in part:
                        if pnt:
                            if pntCount > 0:
                                cumLength += Dist(startx, starty, pnt.X, pnt.Y) 
                            PointX.append(pnt.X)
                            PointY.append(pnt.Y)
                            PointZ.append(pnt.Z)
                            LengthfromStart.append(cumLength)

                            startx = pnt.X
                            starty = pnt.Y
                            pntCount += 1

                ##Step 1: Split the cross section by the lowest points
                pointZArr = (np.array(PointZ)*100).astype(int)
                min_Z = min(pointZArr)

                ##Seperate the pointX to two arrays based on the lowest elevation
                array = np.append(pointZArr, np.inf)  # padding so we don't lose last element
                pointXarr = np.append(PointX, np.inf)  # padding so we don't lose last element
                pointYarr = np.append(PointY, np.inf)  # padding so we don't lose last element
                floatZarr = np.append(PointZ, np.inf)  # padding so we don't lose last element
                LengthArr = np.append(LengthfromStart, np.inf)
                

                split_indices = np.where(array == min_Z)[0]
                splitarray = np.split(array, split_indices + 1)
                splitpointXarr = np.split(pointXarr, split_indices + 1)
                splitpointYarr = np.split(pointYarr, split_indices + 1)
                splitpointZarr = np.split(floatZarr, split_indices + 1)
                splitlengtharr = np.split(LengthArr, split_indices + 1)

                ##Cut the cross section by the lowest point and then to cut the highest points into each half
                k = 0
                for subarray in splitarray:
                    #arcpy.AddMessage("side: #" + str(k))
                    if len(subarray) > 5: ##the half profile should at least have 5 points
                        half_pointZarr = subarray[:-1]
                        #arcpy.AddMessage(half_pointZarr/100)

                        subpointXarr = splitpointXarr[k]
                        half_pointXarr = subpointXarr[:-1]
                        
                        subpointYarr = splitpointYarr[k]
                        half_pointYarr = subpointYarr[:-1]

                        sublengtharr = splitlengtharr[k]
                        half_lengtharr = sublengtharr[:-1]

                        z_max = max(half_pointZarr)
                        z_min = min(half_pointZarr)
                        #arcpy.AddMessage(z_min)
                        idx = np.where(half_pointZarr == z_max)[0][0]

                        ##Record the X and Y coordinates for the highest points
                        if (k ==0 or k==len(splitarray)-1):
                            X_coord.append(half_pointXarr[idx])
                            Y_coord.append(half_pointYarr[idx])
                            pntType.append(1)  ##1: highest points
                            FID.append(oid)
                            height = int((z_max - z_min)/100 + 0.5)
                            #arcpy.AddMessage(height)
                            Height.append(height)
                            side.append(k)
                            width = int(abs(half_lengtharr[-1] - half_lengtharr[0])+ 0.5)
                            Length.append (width)

                        if "convex" in AdjustProfile: 
                            #arcpy.AddMessage("Cut the cross sections by the largest convex points on each side...")
                            #init_height = int(min_height)
                            if (half_pointZarr[0] > half_pointZarr[-1]): ##Left side of the profile
                                validpointZarr = half_pointZarr[idx:]
                                validpointXarr = half_pointXarr[idx:]
                                validpointYarr = half_pointYarr[idx:]
                                validlengtharr = half_lengtharr[idx:]
                                validlengtharr = validlengtharr - min(validlengtharr) ##normalize the length values
                            else: ##Right-side of the profile; reverse the order of the array
                                validpointZarr = np.flip(half_pointZarr[:idx+1])
                                validpointXarr = np.flip(half_pointXarr[:idx+1])
                                validpointYarr = np.flip(half_pointYarr[:idx+1])
                                validlengtharr = np.flip(half_lengtharr[:idx+1])
                                validlengtharr = max(validlengtharr) - validlengtharr ##normalize the length values

                            idx = 0
                            dist = 1000

                            looplengtharr = validlengtharr[idx:]
                            looppointZarr = validpointZarr[idx:]
                            looppointXarr = validpointXarr[idx:]
                            looppointYarr = validpointYarr[idx:]
                            loop_height = 10*min_height

                            while (True):
                                max_idx, max_dist  = turning_points_RDP(looplengtharr, looppointZarr/100) ###, 1, int(cellsize_float)*3) ##only top 1 turing points should be enough
                                idx = max_idx
                                dist = max_dist
                                loop_height = (max(looppointZarr[idx:]) - z_min)/100
                                width = looplengtharr[-1] - looplengtharr[idx]

                                adj_min_height = max(width/min_width * min_height, min_height) ##Adjust the min_height based on the width

                                if (dist > 20 and loop_height > adj_min_height and idx > 0):
                                    X_coord.append(looppointXarr[idx])
                                    Y_coord.append(looppointYarr[idx])
                                    pntType.append(2)  ##1: convex points
                                    FID.append(oid)
                                    Height.append(int(loop_height+0.5))
                                    side.append(k)
                                    Length.append (int(width+0.5))

                                    looplengtharr = looplengtharr[idx:]
                                    looppointZarr = looppointZarr[idx:]
                                    looppointXarr = looppointXarr[idx:]
                                    looppointYarr = looppointYarr[idx:]
                                else:
                                    break
                    k += 1        

        bnd_points = arcpy.CreateFeatureclass_management(temp_workspace, "bnd_points","POINT", "","","", spatialref)
        arcpy.AddField_management(bnd_points, 'PntID', 'Long', 6) 
        arcpy.AddField_management(bnd_points, 'SectionID', 'Long', 6) 
        arcpy.AddField_management(bnd_points, 'Side', 'Long', 6) 
        arcpy.AddField_management(bnd_points, 'PntType', 'String', 10) 
        arcpy.AddField_management(bnd_points, 'Height', 'Long', 10) 
        arcpy.AddField_management(bnd_points, 'Width', 'Long', 10) 

        new_point_cursor = arcpy.da.InsertCursor(bnd_points, ('SHAPE@', 'PntID', 'SectionID', 'Side','PntType', 'Height', 'Width'))
        for i in range(len(X_coord)):
            pnt = arcpy.Point(X_coord[i],Y_coord[i])
            ptype = pntType[i]
            if ptype == 1:
                ptypeStr = "Highest"
            else:
                ptypeStr = "Convex"
            new_point_cursor.insertRow([pnt, i, FID[i], side[i], ptypeStr, Height[i], Length[i]])
        del new_point_cursor

        if OutputConvexPoints != "":
            arcpy.CopyFeatures_management(bnd_points, OutputConvexPoints)            
        

        arcpy.management.SplitLineAtPoint(singlepartlines, bnd_points, temp_workspace + "\\split_profiles", "1 Meters")
        fieldmappings = arcpy.FieldMappings()
        fieldmappings.addTable(singlepartlines)
        ##Should use the lowest points for the spatial join
        arcpy.SpatialJoin_analysis(temp_workspace + "\\split_profiles", flowline3dpoints, final_cross_sections, "JOIN_ONE_TO_ONE", "KEEP_COMMON", fieldmappings, "INTERSECT", "1 Meters", "#")
    else:
        arcpy.CopyFeatures_management(singlepartlines, final_cross_sections) 

    if b_divide:
        arcpy.AddMessage("Step 5: Divide cross sections by the streamlines...")
        splitted_perps = temp_workspace + "\\splitted_perps"
        arcpy.SplitLineAtPoint_management(final_cross_sections, flowline3dpoints, splitted_perps, '1 Meters')
        ##Delete the cross section line that is less than 90% of the width
        with arcpy.da.UpdateCursor(splitted_perps, 'SHAPE@LENGTH') as cursor:
            for row in cursor:
                if row[0] < min_width:
                    cursor.deleteRow()
        del row, cursor
        arcpy.CopyFeatures_management(splitted_perps, out_cross_sections)
    else:
        with arcpy.da.UpdateCursor(final_cross_sections, 'SHAPE@LENGTH') as cursor:
            for row in cursor:
                if row[0] < min_width:
                    cursor.deleteRow()
        del row, cursor
        arcpy.CopyFeatures_management(final_cross_sections, out_cross_sections)

####-------Start the main program-----------------------####
if __name__ == '__main__':

    #Define input data this is the core data
    BedDEM = arcpy.GetParameterAsText(0)
    inputflowline = arcpy.GetParameterAsText(1)
    constrainboundary = arcpy.GetParameterAsText(2)
    eraseAreas = arcpy.GetParameterAsText(3)
    spacing=int(arcpy.GetParameter(4))
    half_width = int(arcpy.GetParameter(5))

    AdjustProfile = arcpy.GetParameterAsText(6)
    min_width = int(arcpy.GetParameter(7))
    min_height = arcpy.GetParameter(8)
    b_divide = arcpy.GetParameter(9)
    out_cross_sections=arcpy.GetParameterAsText(10)
    OutputConvexPoints  = arcpy.GetParameterAsText(11) ##Input turning points or cross sections around the outlet points
    OutputFolder = arcpy.GetParameterAsText(12)


    arcpy.Delete_management(temp_workspace)

    singlepartlines = CreateCrossSections(BedDEM, inputflowline, constrainboundary, eraseAreas, spacing, half_width, AdjustProfile, min_width, min_height, b_divide, out_cross_sections, OutputConvexPoints)

    if OutputFolder != "":
        arcpy.AddMessage("Step 6: Save cross-sectional plots...")
        arcpy.AddField_management(out_cross_sections, "ProfilePlot", "TEXT", 20)
        arcpy.AddField_management(out_cross_sections, "ProfileID", "Long", 10)
        arcpy.CalculateField_management(out_cross_sections,"ProfileID",str("!FlowPntID!"),"PYTHON_9.3")

        profile3D = temp_workspace + "\\profile3D"
        arcpy.InterpolateShape_3d(BedDEM, out_cross_sections, profile3D)

        plot_list = []
        FID_list = []
        with arcpy.da.SearchCursor(profile3D, ["ProfileID", "SHAPE@", "SHAPE@LENGTH"]) as cursor:
            for row in cursor: ##Loop for each line
                LengthfromStart = []
                PointZ = []
                fcID = row[0]
                FID_list.append(fcID)
                
                for part in row[1]:
                    cumLength = 0
                    pntCount = 0
                    for pnt in part:
                        if pnt:
                            if pntCount > 0:
                                cumLength += Dist(startx, starty, pnt.X, pnt.Y) 
                            PointZ.append(pnt.Z)
                            LengthfromStart.append(cumLength)

                            startx = pnt.X
                            starty = pnt.Y
                            pntCount += 1

                #Save the topographic plots
                fig, ax = plt.subplots()
                ax.plot(LengthfromStart, PointZ)
                ax.set_title(f'Cross Section: # ProfileID: {fcID}')
                ax.set_xlabel('Distance (m)')
                ax.set_ylabel('Elevation (m)')
                filename = OutputFolder + "\\ProfileID_" + str(fcID)+".png"
                fig.savefig(filename, dpi=300, bbox_inches='tight')
                plt.close(fig)  # Close the figure to save computer processing
                plotlink = "file:///" + filename
                plot_list.append(plotlink)

        del row, cursor

        fields = ("ProfileID", "ProfilePlot")
        with arcpy.da.UpdateCursor(out_cross_sections, fields) as cursor:
            for row in cursor:
                fid = FID_list.index(row[0])
                row[1] = plot_list[fid]
                
                #update cursor
                cursor.updateRow(row)

        del row, cursor

    arcpy.Delete_management(temp_workspace) 
