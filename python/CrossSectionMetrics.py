#-------------------------------------------------------------------------------
# Name: CrossSectionMetrics.py
#
# Purpose:
# This tool derives metrics for cross sections (two-sided profiles). The inputs include
# the DEM, cross sections (can be a file or digitized on screen), the method to refine
# cross sections, the minimum cross section width between two sides, and the minimum height
# above the valley floor when refining cross sections. The outputs include the refined cross
# sections with derived metrics and optional points to refine the cross sections and half
# cross-section metrics. Because users may provide their own digitized cross sections, this
# tool also provides the same three options described in the previous tool to refine the cross
# sections. To help determine the appropriate cross sections, the tool also generates the highest
# and convex points on both sides of the valley as an optional output. To help facilitate
# comparisons of morphologies on each side of the valley, the tool also provides an option
# to save the half cross-section metrics. Similar to the previous tool, this tool also provides
# an option to specify a folder in which to save the plots of cross-sectional profiles for
# further analysis.
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
import math
import time
import numpy as np
from scipy.optimize import curve_fit
from scipy import optimize
import matplotlib.pyplot as plt

arcpy.env.overwriteOutput = True
arcpy.env.XYTolerance= "0.01 Meters"

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

# Polynomial Regression
def polyfit(x, y, degree):
    results = {}

    coeffs = np.polyfit(x, y, degree)

     # Polynomial Coefficients
    results['polynomial'] = coeffs.tolist()

    # r-squared
    p = np.poly1d(coeffs)
    # fit values, and mean
    yhat = p(x)                         # or [p(z) for z in x]
    ybar = np.sum(y)/len(y)          # or sum(y)/len(y)
    ssreg = np.sum((yhat-ybar)**2)   # or sum([ (yihat - ybar)**2 for yihat in yhat])
    sstot = np.sum((y - ybar)**2)    # or sum([ (yi - ybar)**2 for yi in y])
    results['determination'] = ssreg / sstot

    return results

def k_curve (x, c):
    return (1-x) * np.exp(c * x)

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

#---------------------------------------------------------------------------------------------------------------
# This function calculates the distance between two points
#--------------------------------------------------------------------------------------------------------------- 
def Dist(x1,y1,x2,y2):
    return math.sqrt(math.pow(math.fabs(x1-x2),2)+math.pow(math.fabs(y1-y2),2))

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
        arcpy.FlipLine_edit("lyrLines")  ##Need to change to lyrLines
        arcpy.SelectLayerByAttribute_management("lyrLines", "CLEAR_SELECTION")

    arcpy.DeleteField_management (line, "Flip")
    arcpy.Delete_management(line3d)

###rdp only positive distance!!! for turning point detection
def Knickpoints_rdp(points, epsilon, turn_points, dists):
    # get the start and end points
    start = np.tile(np.expand_dims(points[0], axis=0), (points.shape[0], 1))
    end = np.tile(np.expand_dims(points[-1], axis=0), (points.shape[0], 1))
    linedist = Dist(start[0][0],start[0][1],end[0][0],end[0][1])

    dist_point_to_line = np.cross(end - start, points - start, axis=-1) / np.linalg.norm(end - start, axis=-1)
    arcpy.AddMessage("Distance to line")
    arcpy.AddMessage(dist_point_to_line)

    max_idx = np.argmax(dist_point_to_line)
    arcpy.AddMessage(max_idx)
    max_value = dist_point_to_line[max_idx]##/linedist
    arcpy.AddMessage(max_value)
    max_value2 = max(dist_point_to_line)
    arcpy.AddMessage(max_value2)

    if abs(max_value) > epsilon:  ##the distance is at least 1 m from the line
        if max_value > 0:
             turn_points.append(points[max_idx])
             dists.append(max_value)

        partial_results_left = Knickpoints_rdp(points[:max_idx+1], epsilon, turn_points,dists)
        partial_results_right = Knickpoints_rdp(points[max_idx:], epsilon, turn_points,dists)


def turning_points(array):
    ''' turning_points(array) -> min_indices, max_indices
    Finds the turning points within an 1D array and returns the indices of the minimum and 
    maximum turning points in two separate lists.
    '''
    idx_max, idx_min = [], []
    if (len(array) < 3): 
        return idx_min, idx_max

    NEUTRAL, RISING, FALLING = range(3)
    def get_state(a, b):
        if a < b: return RISING
        if a > b: return FALLING
        return NEUTRAL

    ps = get_state(array[0], array[1])
    begin = 1
    for i in range(2, len(array)):
        s = get_state(array[i - 1], array[i])
        if s != NEUTRAL:
            if ps != NEUTRAL and ps != s:
                if s == FALLING: 
                    idx_max.append((begin + i - 1) // 2)
                else:
                    idx_min.append((begin + i - 1) // 2)
            begin = i
            ps = s
    return idx_min, idx_max


def turning_points_RDP(streamLength, streamZ): ###, turning_points = 10, cluster_radius = 200):
    ##use power law to fit the cross section and then calculate the different between Z and Z (curve)

    #LengthfromStart.reverse()
    max_Z = max(streamZ)
    min_Z = min(streamZ)
    max_len = streamLength[-1]
    min_len = streamLength[0]
    normalH = np.array([(y - min_Z)/(max_Z - min_Z) for y in streamZ]) ##+ 0.01
    normalLen = np.array([(y-min_len) /(max_len-min_len) for y in streamLength]) ##+ 0.01

    try: 
        fit_results = k_curve_fit(normalLen, normalH)
        c = fit_results[0]
        R2 = fit_results[1]
        residual = fit_results[2]
    except:
        c = 100

    if c < 0:
        max_idx = np.argmax(residual)
        max_value = residual[max_idx] * (max_Z - min_Z) ##/linedist
    else:
        points = np.concatenate([np.expand_dims(streamLength, axis=-1), np.expand_dims(streamZ, axis=-1)], axis=-1)
        start = np.tile(np.expand_dims(points[0], axis=0), (points.shape[0], 1))
        end = np.tile(np.expand_dims(points[-1], axis=0), (points.shape[0], 1))
        dist_point_to_line = np.cross(end - start, points - start, axis=-1) / np.linalg.norm(end - start, axis=-1)

        max_idx = np.argmax(dist_point_to_line)
        max_value = dist_point_to_line[max_idx]##/linedist

    return max_idx, max_value

##Main program
# Script arguments
InputDEM = arcpy.GetParameterAsText(0)
InputProfiles = arcpy.GetParameterAsText(1)
AdjustProfile = arcpy.GetParameterAsText(2)
min_width = arcpy.GetParameter(3)
min_height = arcpy.GetParameter(4)
OutputProfileMetrics  = arcpy.GetParameterAsText(5) 
OutputConvexPoints  = arcpy.GetParameterAsText(6)
OutputHalfProfileMetrics  = arcpy.GetParameterAsText(7) 
OutputFolder = arcpy.GetParameterAsText(8)
#environments

spatialref=arcpy.Describe(InputProfiles).spatialReference #get spat ref from input
arcpy.env.outputCoordinateSystem = spatialref #output coordinate system is taken from spat ref
arcpy.env.overwriteOutput = True #every new created file with the same name as an already existing file will overwrite the previous file
arcpy.env.XYTolerance= "1 Meters"
arcpy.env.scratchWorkspace=arcpy.env.scratchGDB #define a default folder/database where intermediate product will be stored

cellsize = arcpy.GetRasterProperties_management(InputDEM,"CELLSIZEX")
cellsize_float = float(cellsize.getOutput(0)) # use float cell size


arcpy.Delete_management(temp_workspace) ### Empty the in_memory
profile3D = temp_workspace + "\\profile3D"
arcpy.InterpolateShape_3d(InputDEM, InputProfiles, profile3D)


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

lowest_points = arcpy.CreateFeatureclass_management(temp_workspace, "lowest_points","POINT", "","","", spatialref)
arcpy.AddField_management(lowest_points, 'PntID', 'Long', 6) 

new_point_cursor = arcpy.da.InsertCursor(lowest_points, ('SHAPE@', 'PntID'))
for i in range(len(lowest_X_coord)):
    pnt = arcpy.Point(lowest_X_coord[i],lowest_Y_coord[i])
    new_point_cursor.insertRow([pnt, i])
del new_point_cursor        

if len(AdjustProfile) > 10:  ##Refine profiles 
    ##check the knickpoint idetification tool in AutoCirque
    ##This process will generate a new set of profiles for the following calculations
    #arcpy.InterpolateShape_3d(InputDEM, InputProfiles, temp_workspace + "\\profile3D")
    if "convex" in AdjustProfile:
        arcpy.AddMessage("Cut the cross sections by the largest convex points on each side...")
    else:
        arcpy.AddMessage("Cut the cross sections by the highest points on each side...")
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
                if len(subarray) > 5: ##the half profile should at least have 5 points
                    half_pointZarr = subarray[:-1]

                    subpointXarr = splitpointXarr[k]
                    half_pointXarr = subpointXarr[:-1]
                    
                    subpointYarr = splitpointYarr[k]
                    half_pointYarr = subpointYarr[:-1]

                    sublengtharr = splitlengtharr[k]
                    half_lengtharr = sublengtharr[:-1]

                    z_max = max(half_pointZarr)
                    z_min = min(half_pointZarr)
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

                    #plt.plot( half_lengtharr, half_pointZarr/100)
                    #plt.show()
                    if "convex" in AdjustProfile: 
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

                        #loop = 1
                        while (True):
                            #turning_points_RDP(looplengtharr, looppointZarr/100) ###, 1, int(cellsize_float)*3) ##only top 1 turing points should be enough
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

    #arcpy.CopyFeatures_management(lowest_points, "d:\\temp\\lowest_points.shp") 
    bnd_points = arcpy.CreateFeatureclass_management(temp_workspace + "", "bnd_points","POINT", "","","", spatialref)
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

    ##Check the cutting points and choose the most reasonable coutpoints for each cross section
    SectionArr = np.array(FID)
    HeightArr = np.array(Height)
    WithArr = np.array(Length)
    SideArr = np.array(side)
    pntTypeArr = np.array(pntType)
    X_coordArr = np.array(X_coord)
    Y_coordArr = np.array(Y_coord)
    
    unique_section_arr = np.unique(SectionArr)

    final_FIDs = []
    final_X_coord = []
    final_Y_coord = []
    final_pntType = []
    final_Height = []
    final_side = []
    final_Length = []
    
    for section in unique_section_arr:
        sel_sides = SideArr[SectionArr == section]
        sel_heights = HeightArr[SectionArr == section]
        sel_widths = WithArr[SectionArr == section]
        sel_X_coord = X_coordArr[SectionArr == section]
        sel_Y_coord = Y_coordArr[SectionArr == section]
        sel_pntType = pntTypeArr[SectionArr == section]

        unique_sides = np.unique(sel_sides)
        
        if len(unique_sides) > 1:
            ##Left side
            leftside = unique_sides[0]
            left_heights = sel_heights[sel_sides == leftside]
            left_widths = sel_widths[sel_sides == leftside]
            left_Xs = sel_X_coord[sel_sides == leftside]
            left_Ys = sel_Y_coord[sel_sides == leftside]
            left_pntTypes = sel_pntType[sel_sides == leftside]
            
            rightside = unique_sides[-1]
            right_heights = sel_heights[sel_sides == rightside]
            right_widths = sel_widths[sel_sides == rightside]
            right_Xs = sel_X_coord[sel_sides == rightside]
            right_Ys = sel_Y_coord[sel_sides == rightside]
            right_pntTypes = sel_pntType[sel_sides == rightside]

            min_diff = 100000
            diff_matrix = np.zeros((len(left_heights), len(right_heights)))
            for i in range(len(left_heights)):
                for j in range(len(right_heights)):
                    height_diff = (abs(left_heights[i] - right_heights[j])+1)/(left_heights[i] + right_heights[j]) * 100 ##add 1 to avoid zero difference
                    width_diff = (abs(left_widths[i] - right_widths[j])+1)/(left_widths[i] + right_widths[j]) * 100

                    total_diff = int(height_diff * width_diff * 100)
                    diff_matrix[i][j] = total_diff
                    if total_diff < min_diff:
                        min_diff = total_diff
                    
            x, y = np.where(diff_matrix == min_diff)
            if len(x)>0:
                left_height = left_heights[x[0]]
                left_width = left_widths[x[0]]
                left_X = left_Xs[x[0]]
                left_Y = left_Ys[x[0]]
                left_pntType = left_pntTypes[x[0]]
                
                right_height = right_heights[y[0]]
                right_width = right_widths[y[0]]
                right_X = right_Xs[y[0]]
                right_Y = right_Ys[y[0]]
                right_pntType = right_pntTypes[y[0]]

                ##Add left info            
                final_FIDs.append(section)
                final_side.append(leftside)
                final_X_coord.append(left_X)
                final_Y_coord.append(left_Y)
                final_pntType.append(left_pntType)
                final_Height.append(left_height)
                final_Length.append(left_width)

                ##Add right info
                final_FIDs.append(section)
                final_side.append(rightside)
                final_X_coord.append(right_X)
                final_Y_coord.append(right_Y)
                final_pntType.append(right_pntType)
                final_Height.append(right_height)
                final_Length.append(right_width)


    final_bnd_points = arcpy.CreateFeatureclass_management(temp_workspace + "", "final_bnd_points","POINT", "","","", spatialref)
    arcpy.AddField_management(final_bnd_points, 'PntID', 'Long', 6) 
    arcpy.AddField_management(final_bnd_points, 'SectionID', 'Long', 6) 
    arcpy.AddField_management(final_bnd_points, 'Side', 'Long', 6) 
    arcpy.AddField_management(final_bnd_points, 'PntType', 'String', 10) 
    arcpy.AddField_management(final_bnd_points, 'Height', 'Long', 10) 
    arcpy.AddField_management(final_bnd_points, 'Width', 'Long', 10) 

    new_point_cursor = arcpy.da.InsertCursor(final_bnd_points, ('SHAPE@', 'PntID', 'SectionID', 'Side','PntType', 'Height', 'Width'))
    for i in range(len(final_X_coord)):
        pnt = arcpy.Point(final_X_coord[i],final_Y_coord[i])
        ptype = final_pntType[i]
        if ptype == 1:
            ptypeStr = "Highest"
        else:
            ptypeStr = "Convex"
        new_point_cursor.insertRow([pnt, i, final_FIDs[i], final_side[i], ptypeStr, final_Height[i], final_Length[i]])
    del new_point_cursor

    if OutputConvexPoints != "":
        arcpy.CopyFeatures_management(final_bnd_points, OutputConvexPoints)            
    

    arcpy.management.SplitLineAtPoint(InputProfiles, final_bnd_points, temp_workspace + "\\split_profiles", "1 Meters")
    fieldmappings = arcpy.FieldMappings()
    fieldmappings.addTable(InputProfiles)
    ##Should use the lowest points for the spatial join
    arcpy.SpatialJoin_analysis(temp_workspace + "\\split_profiles", lowest_points, OutputProfileMetrics, "JOIN_ONE_TO_ONE", "KEEP_COMMON", fieldmappings, "INTERSECT", "1 Meters", "#")
else:
    arcpy.CopyFeatures_management(InputProfiles, OutputProfileMetrics) 

if OutputHalfProfileMetrics != "":
    arcpy.management.SplitLineAtPoint(OutputProfileMetrics, lowest_points,OutputHalfProfileMetrics, "1 Meters")
    ##remove the small sections
    min_length = max(cellsize_float * 3, 100)
    with arcpy.da.UpdateCursor(OutputHalfProfileMetrics, 'SHAPE@LENGTH') as cursor:
        for row in cursor:
            if row[0] < min_length:
                cursor.deleteRow()
    del row, cursor

##Derive the whole profile metrics
arcpy.AddMessage("Add profile metric fields...")
Fieldlist=[]
ListFields=arcpy.ListFields(OutputProfileMetrics)

for x in ListFields:
    Fieldlist.append(x.baseName)

if "ProfileID" in Fieldlist:  ## count = 1
    pass
else:
    #add fieds to attribute tables to be populated with calculated values
    arcpy.AddField_management(OutputProfileMetrics, "ProfileID", "LONG", 10)
    arcpy.CalculateField_management(OutputProfileMetrics,"ProfileID",str("!"+str(arcpy.Describe(OutputProfileMetrics).OIDFieldName)+"!"),"PYTHON_9.3")

if OutputFolder != "":
    if "ProfilePlot" in Fieldlist:  ## count = 1
        pass
    else:
        #add fieds to attribute tables to be populated with calculated values
        arcpy.AddField_management(OutputProfileMetrics, "ProfilePlot", "TEXT", 20)

new_fields = ("Length","Height") ## All float variables 4
for field in new_fields:
    if field in Fieldlist:
        pass
    else:
        arcpy.AddField_management(OutputProfileMetrics, field, "DOUBLE", 10, 1)

##Axis variables
new_fields = ("WHRatio", "Asymmetry", "HHRatio", "PI", "V_index") ## All float variables 4
for field in new_fields:
    if field in Fieldlist:
        pass
    else:
        arcpy.AddField_management(OutputProfileMetrics, field, "DOUBLE",10, 3)

  
##axis curve-fit variables 
new_fields = ("Quad_c", "Quad_r2", "VWDR_m", "VWDR_n", "VWDR_r2") ##float variables with high digits
for field in new_fields:
    if field in Fieldlist:
        pass
    else:
        arcpy.AddField_management(OutputProfileMetrics, field, "DOUBLE",10, 4)

arcpy.AddMessage("Derive profile metrics...")
arcpy.InterpolateShape_3d(InputDEM, OutputProfileMetrics, profile3D) 

FID_list = []
PR_list = []
WH_list = []
HH_list = []
asymmetry_list = []
Amp_list = []
length_list = []
v_index_list = []

VWDR_m_list  = []
VWDR_n_list  = []
VWDR_r2_list = []
quad_c_list =  []
quad_r2_list = []

plot_list = []

with arcpy.da.SearchCursor(profile3D, ["ProfileID", "SHAPE@", "SHAPE@LENGTH"]) as cursor:
    i = 0
    line_length = 0
    for row in cursor: ##Loop for each line
        PointX = []
        PointY = []
        LengthfromStart = []
        PointZ = []
        fcID = row[0]
        FID_list.append(fcID)
        cumLength = 0
        line_length = row[2]
        length_list.append(line_length)
        for part in row[1]:
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
        ##Save the cross section plot to outfolder
        if OutputFolder != "":
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

        ##Calculate the PR value  Need to determine the weighted average PR value????
        pointZArr = (np.array(PointZ)*100).astype(int) ##time 100 to make sure that the elevation can be accurate to 0.01 m
        #arcpy.AddMessage(max(LengthfromStart) - line_length)
        min_Z = min(pointZArr)

        array = np.array(pointZArr)  # padding so we don't lose last element
        floatZArr = np.array(PointZ)  # padding so we don't lose last element
        LengthArr = np.array(LengthfromStart)


        split_indices = np.where(array == min_Z)[0]  ##Divide the array into 2 halves based on the minimum Z
        z_min = PointZ[split_indices[0]] ##Get the float z_min
        #length_at_min = PointZ[split_indices[0]] ##Get the float z_min
        splitarray = np.split(floatZArr, split_indices)
        splitLengtharr = np.split(LengthArr, split_indices)

        #total_sections = len(pointZArr) - 1
        max_length = max(LengthfromStart)
        #arcpy.
        #section_len = max(LengthfromStart) / total_sections
        profile_integal = 0
        #WH_ratio = 0
        weights = []
        valley_maxs = []
        valley_heights = []
        v_under_areas = []
        x_under_areas = []
        if len(splitarray) > 1: ##make sure to have at least two sections
            cum_idx = 0
            for idx in range(len(splitarray)):
                cum_idx += len(splitarray[idx])             
                if idx < len(splitarray)-1: ##before the last section
                    half_arr = np.append(splitarray[idx], z_min)
                    half_lengths = np.append(splitLengtharr[idx], LengthArr[cum_idx])
                else:
                    half_arr = splitarray[idx]
                    half_lengths = splitLengtharr[idx]

                z_max = max(half_arr)
                valley_maxs.append(z_max)
                valley_heights.append(z_max - z_min)
                z_mean = sum(half_arr) / len(half_arr)
                pr = (z_mean - z_min) / (z_max - z_min + 0.001) ##to prevent the divide of zero
                weight = (half_lengths[-1]-half_lengths[0])/max_length

                
                weights.append(weight)

                profile_integal += weight * pr

                v_area_under = (z_max - z_min) * (half_lengths[-1]-half_lengths[0]) * 0.5

                x_area_under = (z_mean - z_min) * (half_lengths[-1]-half_lengths[0])
                
                v_under_areas.append(v_area_under)

                x_under_areas.append(x_area_under)
                
        if len(valley_heights) > 1:
            total_area = sum(valley_heights) * max(LengthfromStart) * 0.5
        else:
            total_area = sum(valley_heights) * max(LengthfromStart)

        v_area = total_area - sum(v_under_areas)

        heights = np.array(PointZ) - z_min
        hhratio = min(heights[0], heights[-1])/ max(heights[0], heights[-1])

        HH_list.append(hhratio)
        #arcpy.AddMessage(heights)
        under_area = sum(x_under_areas)

        cross_area =  total_area - under_area
        
        vindex = cross_area / v_area - 1
        #arcpy.AddMessage(("v_index:",vindex))

        v_index_list.append(vindex)
        
        asymmetry  = min(weights[0], weights[-1])/ max(weights[0], weights[-1])
        #weights[0]/sum(weights)
        asymmetry_list.append(asymmetry)
        #arcpy.AddMessage("Total weight is: " + str(total_weight))
        PR_list.append(profile_integal)
        #form_ratio = max(LengthfromStart)/ (max(PointZ) - min(PointZ) + 0.001) ##May need to do the weighted average too!!!
        height = (max(PointZ) - min(PointZ))
        Amp_list.append(height)

        WH_ratio = line_length / height
        WH_list.append(WH_ratio)

        ##Derive VWDR Li et al (2001)
        ##Find the minimum of the Z-max
        max_elev = min(valley_maxs[0], valley_maxs[-1]) ##only consider the leftmost and rightmost sections of the cross section profile
        
        z_min = min(floatZArr)
        if max_elev < (z_min + 10): ##if the valley is only 10 m deep, then use the half profile??
            max_elev = max(valley_maxs)
        ##create the H and W lists
        height_list = []
        WDratio_list = []

        floatZArr2 = np.array(PointZ)
        num = int((max_elev - z_min)/10)  ##use 10m interval for height lists and width list
        for i in range (num):
            elev = min(z_min + (i+1) * 10, max_elev)
            first_index = np.argmax(floatZArr2 < elev)
            last_index = floatZArr2.size - np.argmax(floatZArr2[::-1] < elev) - 1
            ##For the first point X and Y
            pntX1 = PointX[first_index]
            pntY1 = PointY[first_index]
            elev1 = floatZArr2[first_index]
            
            pntX2 = PointX[first_index-1]
            pntY2 = PointY[first_index-1]
            elev2 = floatZArr2[first_index-1]

            pntXstart = pntX1 + (pntX2 - pntX1) / (elev2 - elev1) * (elev - elev1)            
            pntYstart = pntY1 + (pntY2 - pntY1) / (elev2 - elev1) * (elev - elev1)            

            if last_index < len(floatZArr2)-1:
                ##For the last point X and Y
                pntX1 = PointX[last_index]
                pntY1 = PointY[last_index]
                elev1 = floatZArr2[last_index]
                
                pntX2 = PointX[last_index+1]
                pntY2 = PointY[last_index+1]
                elev2 = floatZArr2[last_index+1]

                deltaX = (pntX2 - pntX1) / (elev2 - elev1) * (elev - elev1)
                deltaY = (pntY2 - pntY1) / (elev2 - elev1) * (elev - elev1)
                
                pntXend = pntX1 + (pntX2 - pntX1) / (elev2 - elev1) * (elev - elev1)            
                pntYend = pntY1 + (pntY2 - pntY1) / (elev2 - elev1) * (elev - elev1)            
            else:
                pntXend = PointX[last_index]
                pntYend = PointY[last_index]
            
            width = Dist(pntXstart,pntYstart,pntXend,pntYend)

            height = (elev - z_min)
            wdratio = width / height

            WDratio_list.append(wdratio)
            height_list.append(height)

        ##Derive the power law model fit for the longtitude profile 
        ##Here it is necessary to use only the heights with the values of more than the minimum heights
        HArr = np.array(height_list)
        WDratioArr = np.array(WDratio_list)
        m_list = []
        n_list = []
        R2_list = []
        for ii in range (10): ##upto 100 m cur
            cutoff_height = ii * 20
            validHArr = HArr[HArr > cutoff_height]
            validWDratioArr = WDratioArr[HArr > cutoff_height]

            polyfit_results = polyfit(np.log(np.array(validHArr)), np.log(np.array(validWDratioArr)), 1)
            n = polyfit_results['polynomial'][0]
            m = np.exp(polyfit_results['polynomial'][1])
            R2 = polyfit_results['determination']

            n_list.append(n)
            m_list.append(m)
            R2_list.append(R2)
            
            if (R2 > 0.8) or (len(validHArr) < 5):
                break
            ii += 1

        max_R2 = max(R2_list)
        idx = R2_list.index(max_R2)


        VWDR_m_list.append(m_list[idx])
        VWDR_n_list.append(n_list[idx])
        VWDR_r2_list.append(R2_list[idx])                    

        ##Derive quadratic equation fit for the profile along the width line 03/15/2023
        polyfit_results = polyfit(LengthfromStart,PointZ, 2)
        c = polyfit_results['polynomial'][0] * 100 ##times 100 to enlarge the data
        R2 = polyfit_results['determination']

        quad_c_list.append(c)
        quad_r2_list.append(R2)         

        i += 1

del row, cursor

if OutputFolder != "":
    fields = ("ProfileID", "PI", "WHRatio", "Height", "Quad_c", "Quad_r2", "Asymmetry", "VWDR_m", "VWDR_n", "VWDR_r2", "V_index", "Length", "HHRatio", "ProfilePlot")
else:
    fields = ("ProfileID", "PI", "WHRatio", "Height", "Quad_c", "Quad_r2", "Asymmetry", "VWDR_m", "VWDR_n", "VWDR_r2", "V_index", "Length", "HHRatio")

with arcpy.da.UpdateCursor(OutputProfileMetrics, fields) as cursor:
    for row in cursor:
        try:
            fid = FID_list.index(row[0])
            row[1] = f'{PR_list[fid]:.2f}'
            row[2] = f'{WH_list[fid]:.2f}'
            row[3] = f'{Amp_list[fid]:.1f}'
            row[4] = f'{quad_c_list[fid]:.4f}'
            row[5] = f'{quad_r2_list[fid]:.3f}'
            row[6] = f'{asymmetry_list[fid]:.3f}'
            row[7] = VWDR_m_list[fid]
            row[8] = VWDR_n_list[fid]
            VWDR_r2 = VWDR_r2_list[fid]
            row[9] = f'{VWDR_r2: .3f}'
            row[10] = f'{v_index_list[fid]:.3f}'
            row[11] = f'{length_list[fid]:.1f}'
            row[12] = f'{HH_list[fid]:.3f}'
            if OutputFolder != "":
                row[13] = plot_list[fid]
            
            #update cursor
            cursor.updateRow(row)
        except:
            arcpy.AddMessage("There is an error in the calculation. Move to the next one")
            pass

del row, cursor

##Derive the half valley profile metrics
if OutputHalfProfileMetrics != "":
    #a "list" where the name of fields from the attributed table are copied in
    arcpy.AddMessage("Derive half valley profile metrics...")
    arcpy.AddMessage("Add half profile metric fields...")
    Fieldlist=[]
    ListFields=arcpy.ListFields(OutputHalfProfileMetrics)

    for x in ListFields:
        Fieldlist.append(x.baseName)

    if "ProfileID" in Fieldlist:  ## count = 1
        pass
    else:
        #add fieds to attribute tables to be populated with calculated values
        arcpy.AddField_management(OutputHalfProfileMetrics, "ProfileID", "LONG", 10)
        arcpy.CalculateField_management(OutputHalfProfileMetrics,"ProfileID",str("!"+str(arcpy.Describe(OutputHalfProfileMetrics).OIDFieldName)+"!"),"PYTHON_9.3")
        
    new_fields = ("Length","Height") ## All float variables 4
    for field in new_fields:
        if field in Fieldlist:
            pass
        else:
            arcpy.AddField_management(OutputHalfProfileMetrics, field, "DOUBLE", 10, 1)

    ##Axis variables
    new_fields = ("WHRatio", "Closure", "PI", "SCI", "NCI", "Aspect","Gradient") ## All float variables 4
    for field in new_fields:
        if field in Fieldlist:
            pass
        else:
            arcpy.AddField_management(OutputHalfProfileMetrics, field, "DOUBLE",10, 2)

    ##axis curve-fit variables 
    new_fields = ("Exp_a","Exp_b","Exp_r2","Pow_a", "Pow_b", "Pow_r2","Kcurve_c","Kcurve_r2","SL", "SL_r2") ##float variables with high digits
    for field in new_fields:
        if field in Fieldlist:
            pass
        else:
            arcpy.AddField_management(OutputHalfProfileMetrics, field, "DOUBLE",10, 4)

    ##Check the direction and flip the length from low to high elevations
    arcpy.AddMessage("Check profile direction and flip it from low to high elevations if necessary...")
    Check_If_Flip_Line_Direction(OutputHalfProfileMetrics, InputDEM)

    arcpy.AddMessage("Derive half profile metrics...")
    arcpy.InterpolateShape_3d(InputDEM, OutputHalfProfileMetrics, temp_workspace + "\\profile3D") 

    FID_list = []
    PI_list = []
    HLAsp_list = []
    P_clos_list = []
    Amplitude_list = []
    profgrad_list = []
    length_list = []
    WH_list = []

    exp_a_list = []
    exp_b_list = []
    exp_r2_list = []

    pow_a_list = []
    pow_b_list = []
    pow_r2_list = []

    kcurve_c_list = []
    kcurve_r2_list = []
    SL_list = []
    SL_r2_list = []
    sci_list = []
    nci_list = []

    with arcpy.da.SearchCursor(temp_workspace + "\\profile3D", ["ProfileID", "SHAPE@", "SHAPE@Length"]) as cursor:
        i = 0
        for row in cursor: ##Loop for each line
            #arcpy.AddMessage("Profile #" + str(i+1))
            PointX = []
            PointY = []
            LengthfromStart = []
            PointZ = []
            FID_list.append(row[0])
            lineLength = float(row[2])
            length_list.append(lineLength)
            cumLength = 0
            for part in row[1]:
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
            ##Calculate the HI value
            #arcpy.AddMessage(len(PointZ))
            max_Z = max(PointZ)
            min_Z = min(PointZ)
            mean_Z = sum(PointZ) / len(PointZ)
            PI = (mean_Z - min_Z) / (max_Z - min_Z)+ 0.001 ##add 0.001 to avoid the divide of zero
            PI_list.append(PI)

            height = max_Z - min_Z
            Amplitude_list.append(height)

            whratio = lineLength / height
            WH_list.append(whratio)
            
            gradient = 180.0/math.pi * math.atan((max_Z - min_Z)/max(LengthfromStart))

            profgrad_list.append(gradient)

            ##Calculate the HL-Aspect
            dx  = PointX[0] - PointX[-1]
            dy  = PointY[0] - PointY[-1]

            aspect = 180.0/math.pi * math.atan2(dy, dx)
            if aspect < 90:
                adj_aspect = 90.0 - aspect
            else:
                adj_aspect = 360 + 90.0 - aspect
            HLAsp_list.append(adj_aspect)

            ##Derive the exponential model fit for the longtitude profile 03/15/2023
            pointZArr = np.array(PointZ)

            #HArr = np.array(pointH)
            HArr = pointZArr - min(pointZArr)
            max_H = max(HArr)
            norm_HArr = HArr / max_H

            LenArr = np.array(LengthfromStart)
            max_len = max(LengthfromStart)
            norm_lenArr = LenArr / max_len

            valid_norm_HArr = norm_HArr[np.logical_and(HArr > 0, LenArr > 0)]
            valid_norm_lenArr = norm_lenArr[np.logical_and(HArr > 0, LenArr > 0)]

            num = len(valid_norm_HArr)
            V_HArr = np.zeros(num)
            for ii in range(num):
                V_HArr[ii] = ii/num

            offsets = valid_norm_HArr - V_HArr
            nci = np.median(offsets)
            nci_list.append(nci)
            sci = 1-2*PI
            sci_list.append(sci)
   

            ##Do the normalized regression!!!         
            try:
                polyfit_results = polyfit(valid_norm_lenArr, np.log(valid_norm_HArr), 1)
                b = polyfit_results['polynomial'][0]
                a = np.exp(polyfit_results['polynomial'][1])
                R2 = polyfit_results['determination']
            except:
                #arcpy.AddMessage("There is an error!")
                b = -999
                a = -999
                R2 = -999

            exp_a_list.append(a)
            exp_b_list.append(b)
            exp_r2_list.append(R2)

            ##Derive the power law model fit for the longtitude profile
            ##Do the normalized regression!!! 
            try:
                polyfit_results = polyfit(np.log(valid_norm_lenArr), np.log(valid_norm_HArr), 1)
                b = polyfit_results['polynomial'][0]
                a = np.exp(polyfit_results['polynomial'][1])
                R2 = polyfit_results['determination']
            except:
                b = -999
                a = -999
                R2 = -999

            pow_a_list.append(a)
            pow_b_list.append(b)
            pow_r2_list.append(R2)
          
            ###Calculate the profile closure
            startx = np.array(LengthfromStart[0:-1])
            endx = np.array(LengthfromStart[1:])
            startz = np.array(PointZ[0:-1])
            endz = np.array(PointZ[1:])
            dzdx = (endz - startz)/(endx - startx)

            slopes = 180/np.pi * np.arctan(dzdx)

            #arcpy.AddMessage(slopes)
            if len(slopes) > 3:
                min_slp = np.min(slopes[0:3])
                max_slp = np.max(slopes[-3:])
            else:
                min_slp = np.min(slopes)
                max_slp = np.max(slopes)

            p_close = max_slp - min_slp
            P_clos_list.append(p_close)

            #K-curve-fit
            max_len = max(LengthfromStart)
            PointZ.reverse()
            LengthfromStart.reverse()
            normalH = np.array([(y - min_Z)/(max_Z - min_Z) for y in PointZ])
            normalLen = np.array([(max_len - y) /(max_len) for y in LengthfromStart])
            
            fit_results = k_curve_fit(normalLen, normalH)
            c = fit_results[0]
            R2 = fit_results[1]

            kcurve_c_list.append(c)
            kcurve_r2_list.append(R2)

            ##derive the SL index?????
            pointZArr = np.array(PointZ)

            #HArr = np.array(pointH)
            HArr = pointZArr - min(pointZArr)
            
            LenArr = np.array(LengthfromStart)
            ReverseLengthArr = max(LenArr) - LenArr

            validHArr = HArr[ReverseLengthArr > 0]
            validLenArr = ReverseLengthArr[ReverseLengthArr > 0]

            try:
                polyfit_results = polyfit(np.log(validLenArr), validHArr, 1)
                sl = polyfit_results['polynomial'][0]
                R2 = polyfit_results['determination']
            except:
                sl = -999
                R2 = -999

            SL_list.append(-sl) ##use the positive value
            SL_r2_list.append(R2)
            
            i += 1

    del row, cursor

    fields = ("ProfileID", "Closure", "PI", "Aspect", "Height", "Gradient", "Exp_a","Exp_b","Exp_r2","Pow_a", "Pow_b", "Pow_r2","Kcurve_c","Kcurve_r2","SL", "SL_r2", "Length", "WHRatio", "SCI", "NCI")

    with arcpy.da.UpdateCursor(OutputHalfProfileMetrics, fields) as cursor:
        for row in cursor:
            try:
                fid = FID_list.index(row[0])
                row[1] = f'{P_clos_list[fid]:.1f}'
                row[2] = f'{PI_list[fid]:.3f}'
                row[3] = f'{HLAsp_list[fid]:.1f}'
                row[4] = f'{Amplitude_list[fid]:.1f}'
                row[5] = f'{profgrad_list[fid]:.1f}'

                row[6] = exp_a_list[fid]
                row[7] = exp_b_list[fid]
                row[8] = f'{exp_r2_list[fid]:.3f}'

                row[9] = pow_a_list[fid]
                row[10] = pow_b_list[fid]
                row[11] = f'{pow_r2_list[fid]:.3f}'

                row[12] = f'{kcurve_c_list[fid]:.4f}'
                row[13] = f'{kcurve_r2_list[fid]:.3f}'
                row[14] = f'{SL_list[fid]:.4f}'
                row[15] = f'{SL_r2_list[fid]:.3f}'

                row[16] = f'{length_list[fid]:.1f}'
                row[17] = f'{WH_list[fid]:.2f}'

                row[18] = f'{sci_list[fid]:.3f}'
                row[19] = f'{nci_list[fid]:.3f}'

                #update cursor
                cursor.updateRow(row)
            except:
                pass
    del row, cursor

arcpy.Delete_management(temp_workspace) ### Empty the in_memory












