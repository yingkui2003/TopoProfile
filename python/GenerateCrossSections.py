#-------------------------------------------------------------------------------
# Name: GenerateCrossSections.py
#
# Purpose:
# This tool generates a set of cross sections along streamlines based on a specified spacing
# and maximum width on each side of the streamline. The inputs of this tool include
# the DEM, streamline (can be a file or digitized on the screen), the specified spacing along
# the streamline, the maximum width on each side of the streamline, and the minimum cross-sectional
# width to save. This tool also includes two optional inputs: one is confining boundaries, such as
# watershed boundaries, and the other is areas of exclusion such as lakes or glaciers that the user
# does not want to be crossed by sections. For the streamline network with tributaries, the tool
# generates the sub-watershed corresponding to each tributary streamline to ensure that cross sections
# are generated only within the corresponding sub-watershed of the streamline. In addition, this tool
# provides three options to further constrain cross-section widths: 1) None (no constraint needed);
# 2) Refine the cross sections by the highest point on each side. This removes sections that potentially
# cross over the highest point on each side of the valley; and 3) Terminated by a point of convexity on
# each side. The method to derive convex points along each side of a cross section is based on the
# algorithm to determine the potential cirque threshold points (convex turning points) along streamline
# profiles (Li and Zhao, 2022). Please check Li and Zhao (2022) for details. This algorithm produces a
# list of major convex points on each side of the initial cross section with a fixed width. The cross
# sections are terminated at the lowest convex points on each side of the valley, above a height threshold.
# Specifically, a minimum height above the valley floor corresponding to a specified minimum cross section
# width is required for this method to avoid selection of a convex point close to the valley bottom.
# The main output is the generated cross-sections along the streamlines. The tool also provides options
# to divide cross sections in half by streamlines, save the highest/convex points along the cross sections,
# and export the topographic plot along each cross-section to a specified folder. The end points and
# topographic plots can be used to manually check, revise, and adjust the cross-sections generated
# by this method.
# 
# Author:      Yingkui Li
# Created:     11/07/2024-03/05/2025
# Copyright:   (c) Yingkui Li 2025
# Department of Geography, University of Tennessee
# Knoxville, TN 37996, USA
#-------------------------------------------------------------------------------
#from __future__ import division
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
import pandas as pd
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
def create_cross_sections(flowlinepoints, field, distance, spacing):  
    spatialref=arcpy.Describe(flowlinepoints).spatialReference
    pointarr=arcpy.da.FeatureClassToNumPyArray(flowlinepoints, ('SHAPE@X', 'SHAPE@Y', "LineID"))
    line_ids = np.array([item[2] for item in pointarr])
    unique_line_ids = np.unique(line_ids)
    perpendiculars = arcpy.CreateFeatureclass_management(temp_workspace, "perpendiculars","POLYLINE", "","","", spatialref)
    arcpy.AddField_management(perpendiculars, field, 'Long', 10)
    new_line_cursor = arcpy.da.InsertCursor(perpendiculars, ('SHAPE@', field))
    ##Loop for each streamline
    for line_id in unique_line_ids:
        arr = pointarr[line_ids == line_id]
        ##Loop for each section to create a cross valley profile
        for i in range(len(arr)-1): 
            startx = arr[i][0]
            starty = arr[i][1]
            endx = arr[i+1][0]
            endy = arr[i+1][1]
            centx = (startx + endx) /2
            centy = (starty + endy) /2
            normal = 1
            if starty==endy or startx==endx:
                if starty == endy:
                    y1 = centy + distance
                    y2 = centy - distance
                    x1 = centx
                    x2 = centx
                if startx == endx:
                    y1 = centy
                    y2 = centy 
                    x1 = centx + distance
                    x2 = centx - distance     
            else:
                ##This is revised codes, much simplified!! 11/09/2024
                offsety = distance * (startx - endx) / spacing
                offsetx = distance * (starty - endy) / spacing
                if (abs(offsetx) > distance) or (abs(offsety) > distance):
                    #arcpy.AddMessage((offsetx,offsety))
                    #arcpy.AddMessage(line_id)
                    normal = 0  
                else:
                    y1 = centy + offsety 
                    y2 = centy - offsety
                    x1 = centx - offsetx
                    x2 = centx + offsetx
            if normal > 0:
                array = arcpy.Array([arcpy.Point(x1,y1),arcpy.Point(x2, y2)])
                polyline = arcpy.Polyline(array)
                new_line_cursor.insertRow([polyline, line_id])
                
    del new_line_cursor
    return perpendiculars
#------------------------------------------------------------------------------------------------------------
# This fuction cleanup and clip the cross sections based on the catchment basins.
#------------------------------------------------------------------------------------------------------------
def clean_cross_sections(perpendiculars, flowlines, IDfield, beddem, constrainboundary, eraseAreas):
    #create minimum bounding geometry, convex hull method"
    if constrainboundary != "":
        mbg = arcpy.MinimumBoundingGeometry_management(constrainboundary, temp_workspace + "\\mbg", "ENVELOPE", "ALL", "","NO_MBG_FIELDS")
    else:
        mbg = arcpy.MinimumBoundingGeometry_management(flowlines, temp_workspace + "\\mbg", "ENVELOPE", "ALL", "","NO_MBG_FIELDS")
    ##create a 300 m buffer around the mbg
    buffer_dist = str(half_width) + " Meter"
    mbg_buf = arcpy.Buffer_analysis(mbg, temp_workspace + "\\mbg_buf", buffer_dist, "", "", "ALL")
    extDEM = ExtractByMask(beddem,mbg_buf)
    arcpy.env.extent = extDEM
    arcpy.env.cellSize = extDEM
    arcpy.env.snapRaster = extDEM ##setup snap raster
    
    ##set the parallelProcessingFactor for large DEMs
    nrow = extDEM.height
    ncol = extDEM.width
    oldPPF = arcpy.env.parallelProcessingFactor
    if (nrow > 1500 or ncol > 1500):
        #arcpy.AddMessage("The DEM has " +str(nrow) + " rows and " + str(ncol) + " columns")
        arcpy.env.parallelProcessingFactor = 0 ##use 0 for large rasters    
    fillDEM =Fill(extDEM)  ##Fill the sink first
    fdir = FlowDirection(fillDEM,"NORMAL") ##Flow direction
    facc = FlowAccumulation(fdir) ##Flow accmulation
    ##covert the flowlines to raster
    streamlink = temp_workspace + "\\streamlink"
    arcpy.conversion.FeatureToRaster(flowlines, IDfield, streamlink)
    outWs = Watershed(fdir, streamlink)
    arcpy.RasterToPolygon_conversion(outWs, temp_workspace + "\\catchments", "SIMPLIFY", "VALUE")
    ##Clip the cross section using the watershed boundary
    if constrainboundary != "":
        arcpy.Clip_analysis (temp_workspace + "\\catchments", constrainboundary, temp_workspace + "\\clipcatchments")
        arcpy.CopyFeatures_management(temp_workspace + "\\clipcatchments", temp_workspace + "\\catchments")
    sel_catchment = temp_workspace + "\\sel_catchment"
    sel_perp = temp_workspace + "\\sel_perp"
    clipedlines = temp_workspace + "\\clipedlines"
    
    linearr=arcpy.da.FeatureClassToNumPyArray(flowlines, (IDfield))
    line_ids = np.array([item[0] for item in linearr])
    unique_line_ids = np.unique(line_ids)
    for line_id in unique_line_ids:
        query = "gridcode = " + str(line_id) 
        arcpy.Select_analysis (temp_workspace + "\\catchments", sel_catchment, query)
        query = IDfield + " = " + str(line_id) 
        arcpy.Select_analysis (perpendiculars, sel_perp, query)
        arcpy.Clip_analysis (sel_perp, sel_catchment, temp_workspace + "\\clip_perp")
        if line_id == unique_line_ids[0]:
            arcpy.CopyFeatures_management(temp_workspace + "\\clip_perp", clipedlines)
        else:
            arcpy.Append_management(temp_workspace + "\\clip_perp", clipedlines, "NO_TEST" )
        
    ##Multipart to singleparts
    singlepartline = arcpy.MultipartToSinglepart_management(clipedlines, temp_workspace + "\\singlepartline")
    fieldmappings = arcpy.FieldMappings()
    fieldmappings.addTable(singlepartline)
    singlepartlines = temp_workspace + "\\singlepartlines"
    arcpy.SpatialJoin_analysis(singlepartline, flowlines, singlepartlines, "JOIN_ONE_TO_ONE", "KEEP_COMMON", fieldmappings, "INTERSECT", "#", "#")
    '''
    ##Step 4: Clean perp lines, so that only the part intersect the flowpoints are preserved
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
    ##Reset parallelProcessingFactor to the default
    arcpy.env.parallelProcessingFactor = oldPPF
    
    if eraseAreas != "":
        fieldmappings = arcpy.FieldMappings()
        fieldmappings.addTable(temp_workspace + "\\singlepartlines")
        arcpy.SpatialJoin_analysis(singlepartlines, eraseAreas, temp_workspace + "\\sel_perps", "JOIN_ONE_TO_ONE", "KEEP_COMMON", fieldmappings, "INTERSECT", "1 Meters", "#")
        arcpy.analysis.Erase(singlepartlines, temp_workspace + "\\sel_perps", temp_workspace + "\\final_perps")
        return temp_workspace + "\\final_perps"
    else:
        return singlepartlines
#------------------------------------------------------------------------------------------------------------
# This fuction is the whole process to reconstruct paleoice based on DEM, input flowlines, ice boundary, and default shear stress
#------------------------------------------------------------------------------------------------------------
def CreateCrossSections(BedDEM, flowline, constrainboundary, eraseAreas, spacing, half_width, AdjustProfile, min_width, min_height, b_divide, out_cross_sections, OutputConvexPoints):
    ####Flow direction and accumulation analysis
    cellsize = arcpy.GetRasterProperties_management(BedDEM,"CELLSIZEX")
    cellsize_float = float(cellsize.getOutput(0)) # use float cell size
    spatialref=arcpy.Describe(flowline).spatialReference #get spat ref from input    
    flowlinepoints = temp_workspace + "\\flowlinepoints"
    cross_sections = temp_workspace + "\\cross_sections"
    arcpy.env.extent = BedDEM
    #Step 1: create points along the flowlines
    arcpy.AddMessage("Step 1: Generate cross valley profiles based on spacing and width...")
    inputflowline = temp_workspace + "\\inputflowline"
    arcpy.CopyFeatures_management(flowline, inputflowline)
    
    PointsFile = temp_workspace + "\\PointsFile"
    PointsWithLineID = temp_workspace + "\\PointsWithLineID"
    intersect_points = temp_workspace + "\\intersect_points"
    #lineID = []
    Points = []
    with arcpy.da.SearchCursor(inputflowline, ['SHAPE@']) as cursor:
        i = 0
        for row in cursor:
            #arcpy.AddMessage("----Processing streamline: #" + str(i+1))
            Length = row[0].length
            #try:
            #    rlist = xrange(1, int(Length), int(spacing))
            #except: ##python 3 xrange is replaced by range
            num = int(Length / spacing)
            #rlist = range(1, num)
            
            for j in range(1, num):            
                Points.append(row[0].positionAlongLine(j*spacing))
                #lineID.append (i)
            i += 1
    del cursor, row
    arcpy.CopyFeatures_management(Points, PointsFile)
    '''
    #Always have errors for a large dataset, using spatial join instead   
    arcpy.AddField_management(PointsFile, 'LineID', 'Long', 10) 
 
    with arcpy.da.UpdateCursor(PointsFile, ['LineID']) as cursor:
        i = 0
        for row in cursor:
            row[0] = lineID[i]
            cursor.updateRow(row)
            i+=1
    del row, cursor  
    '''
    try:
        arcpy.AddField_management(inputflowline, 'LineID', 'Long', 10) ##OFID is not FID, but it is related to FID
    except:
        pass
    arcpy.CalculateField_management(inputflowline,"LineID",str("!"+str(arcpy.Describe(inputflowline).OIDFieldName)+"!"),"PYTHON_9.3")
    arcpy.SpatialJoin_analysis(PointsFile, inputflowline, PointsWithLineID, "JOIN_ONE_TO_ONE", "KEEP_COMMON", "#", "INTERSECT", "#", "#")
    cross_sections = create_cross_sections(PointsWithLineID, 'LineID', half_width, spacing)
    ##Clean Corss sections based on drainage basins
    arcpy.AddMessage("Step 2: Constain cross valley profiles within catchments and confined or excluded areas...")
    cleaned_cross_sections = clean_cross_sections(cross_sections, inputflowline, 'LineID', BedDEM, constrainboundary, eraseAreas)
    
    ##Get the intersect points between the cleaned_cross_sections and flowlines
    arcpy.analysis.Intersect([cleaned_cross_sections, inputflowline], intersect_points, "", "", "point")
    ##refine the cross sections
    arcpy.AddMessage("Step 3: Constrain Cross section widths by the highest or convex points...")
    final_cross_sections = temp_workspace + "\\final_cross_sections"
    if len(AdjustProfile) > 10:  ##Refine profiles 
        profile3D = temp_workspace + "\\profile3D"
        arcpy.InterpolateShape_3d(BedDEM, cleaned_cross_sections, profile3D)
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
        if "convex" in AdjustProfile:
            arcpy.AddMessage("----Cut the cross sections by the convex points on each side...")
        else:
            arcpy.AddMessage("----Cut the cross sections by the highest points on each side...")
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
        
        arcpy.management.SplitLineAtPoint(cleaned_cross_sections, bnd_points, temp_workspace + "\\split_profiles", "1 Meters")
        fieldmappings = arcpy.FieldMappings()
        fieldmappings.addTable(cleaned_cross_sections)
        ##Should use the lowest points for the spatial join
        arcpy.SpatialJoin_analysis(temp_workspace + "\\split_profiles", intersect_points, final_cross_sections, "JOIN_ONE_TO_ONE", "KEEP_COMMON", fieldmappings, "INTERSECT", "1 Meters", "#")
    else:
        arcpy.CopyFeatures_management(cleaned_cross_sections, final_cross_sections) 
    if b_divide:
        arcpy.AddMessage("Step 4: Divide cross sections by the streamlines...")
        splitted_perps = temp_workspace + "\\splitted_perps"
        arcpy.SplitLineAtPoint_management(final_cross_sections, flowlinepoints, splitted_perps, '1 Meters')
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
        arcpy.AddMessage("Save cross-sectional plots...")
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
                csvfile = OutputFolder + "\\ProfileID_" + str(fcID)+".csv"
                df = pd.DataFrame({'Length': np.array(LengthfromStart), 'Elevation': np.array(PointZ)})
                df.to_csv(csvfile, index=False)
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
