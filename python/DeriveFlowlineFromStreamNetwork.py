#-------------------------------------------------------------------------------
# Name: DeriveStreamlinesFromStreamNetwotk.py
#
# Purpose:
# This tool derives a streamline or valley bottom line based on a digital elevation model (DEM)
# and a specified lowest valley boundary or cross section, which can be defined as an input
# feature class file (polylines or polygons) or can be digitized on screen. Other required
# inputs include the minimum source area to start a streamline (km2), the minimum total
# contribution area of a stream before joining another stream (km2), and the minimum tributary
# to main valley ratio (Li and Zhao 2022; Li, 2023). The tool also provides three different ways
# to smooth streamlines: 1) fixed smoothing based on a user-specified smoothing distance; 2)
# varied smoothing based on Kienholz et al. (2014) for glacier studies; and 3) no smoothing.
# The output of this tool includes derived streamlines and optional watershed boundaries
# corresponding to the lowest valley locations. This tool is revised from the PalaeoIce
# toolbox (Li, 2023). 
#
# Author: Dr. Yingkui Li
# Created:     11/07/2024-03/05/2025
# Department of Geography, University of Tennessee
# Knoxville, TN 37996
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
# This function calcuates the 2D distance of two points
#------------------------------------------------------------------------------------------------------------
def distance2points(x1, y1, x2, y2):
    import math
    return math.sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1))

#------------------------------------------------------------------------------------------------------------
# This function derives the angle between two lines based on the length of the three corresponding points
#------------------------------------------------------------------------------------------------------------
def angle(a, b, c):
    import math
    try:
        cosc = (a*a + b*b -c*c)/(2*a*b)
        return math.acos(cosc) / 3.14159 * 180
    except:  
        #arcpy.AddMessage("math error ignored!!!")
        return 180

#------------------------------------------------------------------------------------------------------------
# This function insert new features into an existing feature class (polygon, polyline or multipoint) based on
# a NumPy array. This function is from by online free code.
#------------------------------------------------------------------------------------------------------------
def numpy_array_to_features(in_fc, in_array, geom_fields, id_field):

    # Establish minimum number of x,y pairs to create proper geometry
    min_xy_dict = {'Polygon': 3, 'Polyline': 2, 'Multipoint': 1}
    min_xy_pairs = min_xy_dict[arcpy.Describe(in_fc).shapeType]
 
    if isinstance(geom_fields, list) and len(geom_fields) == 1:
        # Can't access a single field via a list later, extract the
        # only value
        geom_fields = geom_fields[0]
 
    with arcpy.da.InsertCursor(in_fc, ['SHAPE@', id_field]) as cursor:
        unique_array = np.unique(in_array[id_field])  # unique ids
 
        # Iterate through unique sets, get array that matches unique
        # value, convert coordinates to a list and insert via cursor.
        for unique_value in unique_array:
            a = in_array[in_array[id_field] == unique_value]
            if len(a) >= min_xy_pairs:  # skip if not enough x,y pairs
                cursor.insertRow([a[geom_fields].tolist(), unique_value])
            else:
                pass  # skip if not enough x,y pairs
    del cursor
    
    return

#------------------------------------------------------------------------------------------------------------
# The function cleans extrlines based on from and to nodes. If only one to node and no corresponding from node, 
# except for the highest facc section, marking for deletion. The same processes are iterated to remove all extra 
# lines. This is more efficient than other clean line methods based on the intersect of the to node points.
#------------------------------------------------------------------------------------------------------------
def cleanextralineswithtopology(inline,outline, field):
    bflag = 1
    while bflag:
        bflag = 0
        lineArray = arcpy.da.FeatureClassToNumPyArray(inline,['OID@','from_node','to_node', field])
        fromnode = np.array([item[1] for item in lineArray])
        tonode = np.array([item[2] for item in lineArray])
        facc = np.array([item[3] for item in lineArray])
        uniquetonode = np.unique(tonode)
        maxfacc = max(facc)

        lineid = [] ##Record the id for deletion
        for i in range(len(lineArray)):
            linetonode = lineArray[i][2]
            nodecount = np.count_nonzero(uniquetonode == linetonode)
            if nodecount == 1 and not (linetonode in fromnode) and lineArray[i][3] < maxfacc: ###only one tonode except for the highest facc section
                #print "mark for deletion"
                lineid.append(lineArray[i][0])
                bflag = 1

        ##Delete the line marked for deletion
        with arcpy.da.UpdateCursor(inline, "OID@") as cursor:
            for row in cursor:
                if int(row[0]) in lineid:
                    cursor.deleteRow()     
        del cursor, row

    arcpy.CopyFeatures_management(inline, outline)
    return outline

#------------------------------------------------------------------------------------------------------------
# This function smooths the streamlines by adjusting the big turns.
#------------------------------------------------------------------------------------------------------------
def streamline_remove_bigturn(streamline, max_angle, cellsize):
    
    arcpy.SimplifyLine_cartography(streamline, temp_workspace + "\\simply_line", 'POINT_REMOVE', (str(cellsize) + ' Meters'))
    arcpy.FeatureVerticesToPoints_management(temp_workspace + "\\simply_line", temp_workspace + "\\streamline_points", 'All')

    ###Create the new line after removing the outlier points
    spatialref=arcpy.Describe(streamline).spatialReference
    new_line = arcpy.CreateFeatureclass_management(temp_workspace, "new_line","POLYLINE", streamline,"","", spatialref)
    arcpy.AddField_management(new_line, "ORIG_FID", "LONG")
    exist_fields = [f.name for f in arcpy.ListFields(streamline)] #List of current field names in outline layer
    #arcpy.AddMessage(exist_fields)
    fields = exist_fields[2:] ##The first two fields are FID and Geometry
    
    linearray = arcpy.da.FeatureClassToNumPyArray(temp_workspace + "\\simply_line", fields)

    pointarray = arcpy.da.FeatureClassToNumPyArray(temp_workspace + "\\streamline_points", ('SHAPE@X', 'SHAPE@Y','ORIG_FID'))
    line_ids = np.array([item[2] for item in pointarray])
    unique_line_ids = np.unique(line_ids)

    for fid in unique_line_ids:
        arr = pointarray[line_ids == fid]
        ##Methd 2: move this point until the angle is larger than the max_angle
        for row in range(len(arr)):
            if row <(len(arr)-1) and row > 0:#if it is not first or last point of all
                x1 = float(arr[row-1][0])
                y1 = float(arr[row-1][1])
                x = float(arr[row][0])
                y = float(arr[row][1])
                x2 = float(arr[row+1][0])
                y2 = float(arr[row+1][1])
                length1 = distance2points(x1, y1, x, y)
                length2 = distance2points(x2, y2, x, y)
                length  = distance2points(x1, y1, x2, y2)
                pntangle = angle(length1, length2, length)
                if pntangle < max_angle:
                    midx = (x1 + x2)/2
                    midy = (y1 + y2)/2
                    for i in range(5):
                        newx = x + (midx - x) * (i+1) / 5
                        newy = y + (midy - y) * (i+1) / 5
                        length1 = distance2points(x1, y1, newx, newy)
                        length2 = distance2points(x2, y2, newx, newy)
                        pntangle = angle(length1, length2, length)
                        if pntangle > max_angle:
                            arr[row][0] = newx
                            arr[row][1] = newy
                            break
     
        numpy_array_to_features(new_line, arr, ['SHAPE@X', 'SHAPE@Y'], 'ORIG_FID')

    ##Assign field to the new_line
    arcpy.DeleteField_management(new_line, 'ORIG_FID')

    with arcpy.da.UpdateCursor(new_line, fields) as cursor:
        i = 0
        for row in cursor:
            for j in range(len(fields)):
                row[j] = linearray[i][j]
            cursor.updateRow(row)
            i += 1
    del cursor, row

    return new_line

#------------------------------------------------------------------------------------------------------------
# This fuction regroups streamlines to individual ValleyID and dissolve the streamline sections from the top 
# to the lowest points or the confluence points of another streamline. The streamline direction has to be from 
# low to high for each streamline.
#------------------------------------------------------------------------------------------------------------
def Merge_and_Add_ValleyID_by_Topology (streamlines, FaccField, ValleyID, MergeID, outstreamline): 

    streamlinecopy = temp_workspace + "\\streamlinecopy"
    arcpy.CopyFeatures_management(streamlines, streamlinecopy)
    arcpy.AddField_management(streamlinecopy, ValleyID, "Long") #Add GID to the streamline layer
    arcpy.AddField_management(streamlinecopy, MergeID, "Long") #Add MergeID to the streamline layer
    #Calculate Field does not work in ArcGIS Pro, using the following code instead
    with arcpy.da.UpdateCursor(streamlinecopy, [ValleyID, MergeID]) as cursor:
        for row in cursor:
            row[0] = -1
            row[1] = -1
            cursor.updateRow(row)
    del row, cursor

    points=[]
    lines=[]
    facc = []
    valleys = []
    valley_ID = []
    mergeid = []
    mergeused = []
    ##First loop to get the startpoints, lines and the ValleyID list
    with arcpy.da.SearchCursor(streamlinecopy, ["SHAPE@", FaccField, ValleyID, MergeID]) as flows:
        for flow in flows:
            points.append(((flow[0].firstPoint).X, (flow[0].firstPoint).Y))
            lines.append(flow[0])
            facc.append(flow[1])
            valley_ID.append(flow[2]) ##default are all -1 for each line
            mergeid.append(flow[3])
            mergeused.append(0)
    del flow, flows
    
    lenarr = np.array(facc)
    ids = lenarr.argsort()[::-1]

    array=np.array(points)
    ##Second loop: to get the lowest streamline and assign seperated glacier ID
    iceId = 0
    idlist = []
    mergidlist = []
    for i in range(len(array)):
        nTouches=0
        punto=arcpy.Point(array[i][0],array[i][1])
        for a in range (len(array)):
            if a != i: ##don't touch itself
                touches=(punto.touches(lines[a]) or punto.within(lines[a]))
                if touches==True: ##if the start point touches others
                    #arcpy.AddMessage("touched")
                    nTouches = 1
                    break  ##break the loop
        if nTouches==0: ##this is for glaicer iD
            #arcpy.AddMessage("not touched")
            valleys.append(lines[i])
            idlist.append(iceId)
            valley_ID[i] = iceId
            mergeid[i] = iceId
            mergidlist.append(iceId)
            iceId += 1

    ##Loop 3: to test if eachline touches the valleys
    leftover = len(array) - len(valleys)
    while leftover > 0:
        #arcpy.AddMessage("leftover: " + str(leftover))
        for i in range(len(array)):
            nTouches=0
            lineid = ids[i]
            if mergeid[lineid] == -1:
                punto=arcpy.Point(array[lineid][0],array[lineid][1])
                for a in range (len(valleys)): 
                    touches = punto.touches(valleys[a])
                    if touches==True: ##if the start point touches others
                        if valley_ID[lineid] == -1: ##if this line has not been assigned
                            valley_ID[lineid] = idlist[a]
                            valleys.append(lines[lineid])
                            idlist.append(idlist[a])
                            if (mergeused[a] == 0):
                                mergeid[lineid] = mergidlist[a]
                                mergidlist.append(mergidlist[a])
                                mergeused[a] = 1
                            else:
                                ##Start a new mergeID
                                maxmergid = max(mergidlist)
                                mergeid[lineid] = maxmergid + 1
                                mergidlist.append(maxmergid + 1)
                    else:
                        within = punto.within(valleys[a])
                        if within==True: ##if the start point touches others
                            if valley_ID[lineid] == -1: ##if this line has not been assigned
                                valley_ID[lineid] = idlist[a]
                                valleys.append(lines[lineid])
                                idlist.append(idlist[a])
                                ##start a new mergeid with the max mergid + 1
                                maxmergid = max(mergidlist)
                                mergeid[lineid] = maxmergid + 1
                                mergidlist.append(maxmergid + 1)
                        
                ##update leftover
            leftover = len(array) - len(valleys)

    ##Finally add ValleyID to the outstreamline
    i = 0
    with arcpy.da.UpdateCursor(streamlinecopy, [ValleyID, MergeID]) as cursor:
        for row in cursor:
            row[0] = valley_ID[i]
            row[1] = mergeid[i]
            cursor.updateRow(row)
            i += 1
    del row, cursor

    ##Disslove the streamline
    field_treatment = FaccField +" SUM; " + ValleyID +" FIRST"
    arcpy.Dissolve_management(streamlinecopy, outstreamline, MergeID, field_treatment, 'SINGLE_PART', 'UNSPLIT_LINES')

    exist_fields = [f.name for f in arcpy.ListFields(outstreamline)] #List of current field names in outline layer
    for field in exist_fields:
        if field[0:6] == "FIRST_":
            FirstGID = field
        elif field[0:4] == "SUM_":
            MaxFcc = field
    
    new_fields = [FaccField, ValleyID]
    for field in new_fields:
        if field not in exist_fields:
            arcpy.AddField_management(outstreamline, field, "LONG")
    
    fields = [FaccField, ValleyID, MaxFcc, FirstGID]
    with arcpy.da.UpdateCursor(outstreamline, fields) as cursor:
        for row in cursor:
            row[0] = row[2]
            row[1] = row[3]
            cursor.updateRow(row)
    del row, cursor

    arcpy.DeleteField_management(outstreamline,[FirstGID, MaxFcc, MergeID])

    return outstreamline

#------------------------------------------------------------------------------------------------------------
# This function check each line in the line feature and make sure the line is from low elevation to high
# elevation (Glacier streamline needs from low to high elevation in order to reconstruct the paleo ice thickness).
# It is revised from the codes by Pellitero et al.(2016) in GlaRe.
#------------------------------------------------------------------------------------------------------------
def Check_If_Flip_Line_Direction(line, dem):
    cellsize = arcpy.GetRasterProperties_management(dem,"CELLSIZEX")
    cellsize_int = int(float(cellsize.getOutput(0)))
    #arcpy.AddMessage("cellsize_int: " + str(cellsize_int))

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
        arcpy.AddMessage("The number of fliped lines is: " + str(sum(flip_list)))
        arcpy.FlipLine_edit("lyrLines")  ##Need to change to lyrLines
        arcpy.SelectLayerByAttribute_management("lyrLines", "CLEAR_SELECTION")

    arcpy.DeleteField_management (line, "Flip")
    arcpy.Delete_management(line3d)

                    
#------------------------------------------------------------------------------------------------------------
# This fuction smooths the streamlines based on the size of the watershed (flow accumulation)
#------------------------------------------------------------------------------------------------------------
def lineSmooth(inline, outline, smoothfield, cellsize):
    LineID = arcpy.Describe(inline).OIDFieldName
    countResult = arcpy.GetCount_management(inline)
    count = int(countResult.getOutput(0))
    cellarea = float(cellsize) * float(cellsize)
    for i in range (count):
        query = LineID +" = "+ str(i+1)  ##THe FID starts from 1 for the in_memory data; Which is different from shape (FID from 0) and feature class in the geodatabase (objectID from 1)
        #arcpy.AddMessage(query)
        arcpy.Select_analysis(inline, temp_workspace + "\\linesection", query)
        #Determine the smooth tolerance
        areacount = 0
        with arcpy.da.SearchCursor(temp_workspace + "\\linesection", smoothfield) as cursor:
            for row in cursor:
                areacount = row[0]
        del cursor, row

        tolerance = int(areacount * cellarea * 2 / 1.0e6) + 200 ##The function provided by Kienholz et al. (2014) and  James & Carrivick (2016)
        if tolerance > 1000:
            tolerance = 1000
        #arcpy.AddMessage(str(tolerance))
        arcpy.cartography.SmoothLine(temp_workspace + "\\linesection", temp_workspace + "\\linesectsmooth", "PAEK", tolerance)

        if i == 0: ##The first loop
            arcpy.CopyFeatures_management(temp_workspace + "\\linesectsmooth", outline)
        else:
            arcpy.Append_management(temp_workspace + "\\linesectsmooth", outline, "NO_TEST")    
    return outline

#------------------------------------------------------------------------------------------------------------
# This fuction smooths the streamlines based on the size of the watershed (flow accumulation)
#------------------------------------------------------------------------------------------------------------
def lineSmoothFixDistance(inline, outline, smooth_dist):
    arcpy.cartography.SmoothLine(inline, outline, "PAEK", smooth_dist)
    return outline

#------------------------------------------------------------------------------------------------------------
# This fuction is the main program to derive streamlines from stream network.
#------------------------------------------------------------------------------------------------------------
def streamline_from_Stream_Network (InputDEM, InputValleyorCrossSection, StreamThresholdKM2, TributaryThresholdKM2, TributaryRatio, smooth_method, smooth_dis, StreamLine, outWatershed):

    ValleyID = "ValleyID" ##Add a ValleyID for each moriane or cross section

    cellsize = arcpy.GetRasterProperties_management(InputDEM,"CELLSIZEX")
    cellsize_int = int(float(cellsize.getOutput(0)))
    arcpy.env.snapRaster = InputDEM

    StreamThreshold = int(float(StreamThresholdKM2) * 1e6 / (cellsize_int * cellsize_int))
    TributaryThreshold = int(float(TributaryThresholdKM2) * 1e6 / (cellsize_int * cellsize_int))

    ###Step 1: Stream network
    arcpy.AddMessage("Step 1: Stream extraction...")
 
    #Calculate Flowdirection
    fillDEM =Fill(InputDEM)  ##Fill the sink first
    fdir = FlowDirection(fillDEM,"NORMAL") ##Flow direction

    #Calculate Flowaccumulation
    facc = FlowAccumulation(fdir) ##Flow accmulation

    TmpStream = temp_workspace + "\\TmpStream"
    valleyselected = temp_workspace + "\\valleyselected"  ##Set a in_memory file for each moraine feature
    MaxFccTable = temp_workspace + "\\MaxFccTable"
    CleanStream = temp_workspace + "\\CleanStream"
    tmpoutStream = temp_workspace + "\\tmpoutStream"
    smoothline = temp_workspace + "\\smoothline"
    tmpws = temp_workspace + "\\tmpws"
    tmpbuf = temp_workspace + "\\tmpbuf"
    intersect_points = temp_workspace + "\\intersect_points"

    outGreaterThan = Con(facc > StreamThreshold, 1,0)  ##Determine the highest flowaccumuation part
    outStreamLink = StreamLink(outGreaterThan, fdir)
    
    # Process: Stream to Feature
    StreamToFeature(outStreamLink, fdir, TmpStream, "SIMPLIFY")

    FcID = arcpy.Describe(InputValleyorCrossSection).OIDFieldName

    arr=arcpy.da.FeatureClassToNumPyArray(InputValleyorCrossSection, FcID)
    FIds = np.array([item[0] for item in arr])
    count = len(FIds)

    if count < 1:
        arcpy.AddMessage("There is no features in valleys or cross sections! Quit the program!")
        sys.exit()

    bPolyline = True
    fc_type = arcpy.Describe(InputValleyorCrossSection).shapeType
    if fc_type == "Polygon": ##quit if not polyline features
        arcpy.AddMessage("Derive the streamlines within the polygons")
        bPolyline = False


    outstreamline = arcpy.CreateFeatureclass_management(temp_workspace, "outstreamline","POLYLINE","","","",InputValleyorCrossSection)
    arcpy.AddField_management(outstreamline, "Max_Max", "Long") 

    for ivalley in range (count):
        arcpy.AddMessage("Generating streamline(s) for valley #"+str(ivalley + 1)+" of "+str(count) + " valley(s)")

        query = FcID +" = "+str(FIds[ivalley])
        arcpy.Select_analysis(InputValleyorCrossSection, valleyselected, query)

        if bPolyline:
            ##make a small buffer of the cross section to make sure the cross section get the highest fcc
            arcpy.Buffer_analysis(valleyselected, tmpbuf, (str(cellsize_int)+ " Meter"))
            bufID = arcpy.Describe(tmpbuf).OIDFieldName
            
            outZonalStatistics = ZonalStatistics(tmpbuf, bufID, facc, "MAXIMUM") #Find the maximum flowaccumulation point on the moriane feature
        else: ## for polygon input
            outZonalStatistics = ZonalStatistics(valleyselected, FcID, facc, "MAXIMUM")
            
        OutHighestFcc = Con(facc == outZonalStatistics,facc)  ##Determine the highest flowaccumuation part
        
        outSnapPour = SnapPourPoint(OutHighestFcc, facc, 0) ## Just create a pourpoint raster with the same extent of the input DEM
        
        #Calculate Watershed
        outWs = Watershed(fdir, outSnapPour)
        ConOutWs = Con(outWs >= 0, 1)  
        ##Boundary clean
        OutBndCln = BoundaryClean(ConOutWs)

        arcpy.RasterToPolygon_conversion(OutBndCln, tmpws, "NO_SIMPLIFY", "VALUE")

        if bPolyline:
            ##erase the moraines by watershed to see if there is some outside of the watershed
            arcpy.Erase_analysis(valleyselected, tmpws, temp_workspace + "\\valley_outside")
            outsideArr = arcpy.da.FeatureClassToNumPyArray( temp_workspace + "\\valley_outside", 'OID@')

            if len(outsideArr) > 0:
                ws_line = temp_workspace + "\\ws_line"
                arcpy.PolygonToLine_management(tmpws, ws_line)
                arcpy.Append_management(temp_workspace + "\\valley_outside", ws_line, "NO_TEST")                    
                arcpy.FeatureToPolygon_management(ws_line, temp_workspace + "\\wsline_poly", "5 Meters", "NO_ATTRIBUTES")
                arcpy.Dissolve_management(temp_workspace + "\\wsline_poly", tmpws, '#', '#', 'SINGLE_PART')
            
        ##delete potential small polygons
        poly_areaArr = arcpy.da.FeatureClassToNumPyArray(tmpws, 'SHAPE@AREA')
        poly_areas = np.array([item[0] for item in poly_areaArr])
        max_area = np.max(poly_areas)
        if len(poly_areas) > 1:
            with arcpy.da.UpdateCursor(tmpws, 'SHAPE@AREA') as cursor:
                for row in cursor:
                    if int(row[0]) < (max_area - 0.5):
                        cursor.deleteRow()     
            del cursor, row				

        #Get the watershed if required
        if outWatershed !="":
            if ivalley < 1: ##The first loop
                arcpy.CopyFeatures_management(tmpws, outWatershed)
            else:
                arcpy.Append_management(tmpws, outWatershed, "NO_TEST")        


        # Process: Extract by Mask
        #try:
        ExtraFcc = ExtractByMask(facc,tmpws)
        
        # Process: Greater Than
        outGreaterThan = Con(ExtraFcc > StreamThreshold, 1,0)  ##Determine the highest flowaccumuation part
        #need to check if outGreaterThan has the 1 values. If not, no stream will be created
        MaxRasterValue = int((arcpy.GetRasterProperties_management(outGreaterThan, "MAXIMUM").getOutput(0)))
        if MaxRasterValue > 0:
            # Process: Stream Link
            outStreamLink = StreamLink(outGreaterThan, fdir)
            
            # Process: Stream to Feature
            StreamToFeature(outStreamLink, fdir, TmpStream, "SIMPLIFY")

            # Process: Zonal Statistics as Table
            ZonalStatisticsAsTable(outStreamLink, "VALUE", ExtraFcc, MaxFccTable, "DATA", "MAXIMUM")

            # Process: Join Field
            arcpy.JoinField_management(TmpStream, "grid_code", MaxFccTable, "Value", "MAX")  ##Join to get the flow accumulation value
            #arcpy.CopyFeatures_management(TmpStream, "c:\\test\\TmpStream02072023.shp")
            ###This TmpStream already have a to_node in the attibute table, so that it can be used to make the decision
            ##the following is the new to remove and unnecessary lines
            lineArray = arcpy.da.FeatureClassToNumPyArray(TmpStream,['OID@','to_node','MAX'])
            tonode = np.array([item[1] for item in lineArray])
            uniquenode = np.unique(tonode)
            lineid = [] ##Record the id for deletion
            #arcpy.AddMessage("Checking tributary threshold...")
            for i in range(len(uniquenode)):
                selArr = lineArray[tonode == uniquenode[i]]
                fcclist = []
                if len(selArr) > 1: ##Sometimes having more than two end points
                    for j in range(len(selArr)):
                        fcclist.append(selArr[j][2])

                    numselected = len(fcclist)
                    while numselected > 1:
                        minfcc = min(fcclist)
                        if minfcc < TributaryThreshold:
                            for j in range(len(selArr)):
                                if selArr[j][2] == minfcc: ##Remove the smaller fcc one
                                    lineid.append(selArr[j][0])
                                    fcclist.pop(j)
                                    selArr = np.delete(selArr, j)##Remove this one and loop to delete others
                                    numselected = len(selArr)
                                    break ##Only remove one each time

                        else: ##quit the loop if all minfcc are larger than the theshold?? Try to remove more based on the ratio
                            break

            for i in range(len(uniquenode)):
                selArr = lineArray[tonode == uniquenode[i]]
                fcclist = []
                if len(selArr) > 1: ##Sometimes having more than two end points
                    for j in range(len(selArr)):
                        fcclist.append(selArr[j][2])

                    sumfcc = sum(fcclist)
                    fccratio = [x / float(sumfcc) for x in fcclist]

                    for j in range(len(fccratio)):
                        if fccratio[j] < float(TributaryRatio): ##Remove the smaller fcc one
                            lineid.append(selArr[j][0])
            
            
            ##Delete the line marked for deletion
            with arcpy.da.UpdateCursor(TmpStream, "OID@") as cursor:
                for row in cursor:
                    if int(row[0]) in lineid:
                        cursor.deleteRow()     
            del cursor, row				

            ##Clean extralines based on the end points intersection 09/24/2020
            cleanextralineswithtopology(TmpStream,tmpoutStream, 'MAX')  ## clean the extra lines before dissolving

            arcpy.Dissolve_management(tmpoutStream, CleanStream, '#', 'MAX MAX', 'SINGLE_PART', 'UNSPLIT_LINES')

            streamArr = arcpy.da.FeatureClassToNumPyArray(CleanStream, 'OID@')
            if len(streamArr) > 0:
                if "Varied" in smooth_method: 
                    newline = streamline_remove_bigturn(CleanStream, 120, cellsize_int)
                    lineSmooth(newline, smoothline, "Max_Max", cellsize_int)
                    arcpy.Append_management(smoothline, outstreamline, "NO_TEST")
                elif "Fixed" in smooth_method:
                    lineSmoothFixDistance (CleanStream, smoothline, smooth_dis)
                    arcpy.Append_management(smoothline, outstreamline, "NO_TEST")
                else: ##No smooth
                    arcpy.Append_management(CleanStream, outstreamline, "NO_TEST")
            else:
                arcpy.AddMessage("No streamline is created for this feature. It seems that the threshold for a stream is too large!")
        else:
            arcpy.AddMessage("No streamline is created for this feature. It seems that the threshold for a stream is too large!")

    if  bPolyline == False: ##if polygon as the input, clip the streamlines within the polygon
        arcpy.Clip_analysis(outstreamline, InputValleyorCrossSection, temp_workspace + "\\streamline_clip")
        arcpy.CopyFeatures_management(temp_workspace + "\\streamline_clip", outstreamline)        

    ##make sure there are streamlines created
    countResult = arcpy.GetCount_management(outstreamline)
    count = int(countResult.getOutput(0))
    if count < 1:
        arcpy.AddMessage("No streamlines are created for this set of moraine features!!")
        sys.exit()
    ##Merge streamline and add ValleyID
    Check_If_Flip_Line_Direction (outstreamline, fillDEM) ##use fillDEM because the orginal DEM may have problems
    Merge_and_Add_ValleyID_by_Topology (outstreamline, "Max_Max", ValleyID, "MergeID", StreamLine)

##Main program
if __name__ == '__main__':
    # Script arguments
    InputDEM = arcpy.GetParameterAsText(0)
    InputValleyorCrossSection = arcpy.GetParameterAsText(1)
    StreamThresholdKM2 = arcpy.GetParameter(2)
    TributaryThresholdKM2 = arcpy.GetParameter(3)
    TributaryRatio = arcpy.GetParameter(4) ##Define a tributary ratio
    smooth_method = arcpy.GetParameterAsText(5)
    smooth_dis = arcpy.GetParameter(6)
    StreamLine = arcpy.GetParameterAsText(7)
    outWatershed = arcpy.GetParameterAsText(8)

    ##make sure the projection of the glacier outlines is the same with the UTM and the same with the DEM projection 
    spatial_ref_crosssections = arcpy.Describe(InputValleyorCrossSection).spatialReference
    spatial_ref_dem = arcpy.Describe(InputDEM).spatialReference


    if spatial_ref_dem.linearUnitName == "Meter":
        arcpy.AddMessage("The DEM projection is: " + spatial_ref_dem.name)
    else:
        arcpy.AddMessage("The unit of the DEM projection is not in meter. Please re-project the DEM to a projected coordinate system for the analysis!")
        exit()   

    if spatial_ref_crosssections.linearUnitName == "Meter":
        arcpy.AddMessage("The valley cross section projection is: " + spatial_ref_crosssections.name)
    else:
        arcpy.AddMessage("The unit of the valley cross section projection is not in meter. Please re-project it to a projected coordinate system for the analysis!")
        exit()   

    if spatial_ref_dem.name == spatial_ref_crosssections.name:
        arcpy.AddMessage("Both DEM and valley cross section have the same projected coordinate system: " + spatial_ref_dem.name)
    else:
        arcpy.AddMessage("The DEM and valley cross section have different map projections. Please re-project the datasets to the same projection!")
        exit()   
        
    streamline_from_Stream_Network (InputDEM, InputValleyorCrossSection, StreamThresholdKM2, TributaryThresholdKM2, TributaryRatio, smooth_method, smooth_dis, StreamLine, outWatershed)

    ##Delete intermidiate data
    arcpy.Delete_management(temp_workspace) ### Empty the in_memory


   
