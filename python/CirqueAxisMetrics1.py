#-------------------------------------------------------------------------------
# Name:        LongitudinalProfileMetrics
# Purpose:     Derive LongitudinalProfileMetrics based DEM
#
# Author: Yingkui Li
# This program derive cirque related metrics based on cirque outlines and a DEM
# The first step is to determine the cirque threshold points
# The second step is to derive length and width info, as well as the area, and parameters
# The third step is to detive the 3D statistics and hypsometric parameters
# Some of the codes are revised based on the ACME codes by Ramon Pellitero and Matteo Spagnolo 2016
# 
# Created:     05/26/2023
# Copyright:   (c) Yingkui Li 2023
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
            #print "not extension available"
    except:
        raise Exception ("unable to check out extension")
        #print "unable to check out extension"

    try:
        if arcpy.CheckExtension("3D")=="Available":
            arcpy.CheckOutExtension("3D")
        else:
            raise Exception ("not extension available")
            #print "not extension available"
    except:
        raise Exception ("unable to check out extension")
        #print "unable to check out extension"
elif sys.version_info[0] == 3:  ##For ArcGIS Pro
    ArcGISPro = 1
    #pass ##No need to Check
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
    R2 = ssreg / sstot
    if R2 > 1:
        R2 = 1/R2
    return (c, R2)

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
            #arcpy.AddMessage("startZ: " + str(startZ))
            endZ = row[0].lastPoint.Z
            #arcpy.AddMessage("endZ: " + str(endZ))

            if startZ >= endZ:  ##Flip = True use equal in case the start and end point are the same
                flip_list.append(1)
            else:  ##Flip = False
                flip_list.append(0)
            i += 1 

    del cursor
    if i>0:
        del row

    #arcpy.AddMessage(flip_list)
    #arcpy.AddMessage(str(sum(flip_list)))

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


##Main program
# Script arguments
InputDEM = arcpy.GetParameterAsText(0)
InputProfiles = arcpy.GetParameterAsText(1)
b_AdjustProfile = arcpy.GetParameter(2)
#min_height = arcpy.GetParameter(3)
OutputProfileMetrics  = arcpy.GetParameterAsText(3)

#environments
spatialref=arcpy.Describe(InputProfiles).spatialReference #get spat ref from input
arcpy.env.outputCoordinateSystem = spatialref #output coordinate system is taken from spat ref
arcpy.env.overwriteOutput = True #every new created file with the same name as an already existing file will overwrite the previous file
arcpy.env.XYTolerance= "1 Meters"
arcpy.env.scratchWorkspace=arcpy.env.scratchGDB #define a default folder/database where intermediate product will be stored

cellsize = arcpy.GetRasterProperties_management(InputDEM,"CELLSIZEX")
cellsize_float = float(cellsize.getOutput(0)) # use float cell size


arcpy.Delete_management(temp_workspace) ### Empty the in_memory

if b_AdjustProfile: 
    ##Use the highest elevation to cut off the one do not overlap the lowest point
    ##Need to rewrite this part of the program
    
    arcpy.InterpolateShape_3d(InputDEM, InputProfiles, temp_workspace + "\\profile3D")
    try:
        arcpy.AddField_management(temp_workspace + "\\profile3D", "ProfileID", "LONG", 10)
        arcpy.CalculateField_management(temp_workspace + "\\profile3D","ProfileID",str("!"+str(arcpy.Describe(temp_workspace + "\\profile3D").OIDFieldName)+"!"),"PYTHON_9.3")
    except:
        pass
    ##The detailed method below:
    ##Find the highest points and the lowest points for each profile
    ##Split the profile using the hightest points
    ##Spatialjoin with the lowest points to refine the profiles

    arcpy.AddMessage("Refine profiles...")

    max_points = arcpy.CreateFeatureclass_management(temp_workspace + "", "max_points","POINT", "","","", spatialref)
    min_points = arcpy.CreateFeatureclass_management(temp_workspace + "", "min_points","POINT", "","","", spatialref)
    arcpy.AddField_management(max_points, 'PointID', 'Long', 6) 
    arcpy.AddField_management(min_points, 'PointID', 'Long', 6) 

    max_point_cursor = arcpy.da.InsertCursor(max_points, ('SHAPE@', 'PointID'))
    min_point_cursor = arcpy.da.InsertCursor(min_points, ('SHAPE@', 'PointID'))

    max_point_cursor = arcpy.da.InsertCursor(max_points, ('SHAPE@', 'PointID'))
    min_point_cursor = arcpy.da.InsertCursor(min_points, ('SHAPE@', 'PointID'))

    with arcpy.da.SearchCursor(temp_workspace + "\\profile3D", ["ProfileID", "SHAPE@"]) as cursor:
        for row in cursor: ##Loop for each line
            PointX = []
            PointY = []
            PointZ = []
            pntID = row[0]
            for part in row[1]:
                for pnt in part:
                    if pnt:
                        PointX.append(pnt.X)
                        PointY.append(pnt.Y)
                        PointZ.append(pnt.Z)

            pointZArr = np.array(PointZ).astype(int)
            pointXArr = np.array(PointX)
            pointYArr = np.array(PointY)
            max_Z = max(pointZArr)
            min_Z = min(pointZArr)
            maxPointXarr = pointXArr[pointZArr == max_Z]
            maxPointYarr = pointYArr[pointZArr == max_Z]

            for i in range(len(maxPointXarr)):
                pnt = arcpy.Point(maxPointXarr[i],maxPointYarr[i])
                max_point_cursor.insertRow([pnt, pntID])

            minPointXarr = pointXArr[pointZArr == min_Z]
            minPointYarr = pointYArr[pointZArr == min_Z]
            for i in range(len(minPointXarr)):
                pnt = arcpy.Point(minPointXarr[i],minPointYarr[i])
                min_point_cursor.insertRow([pnt, pntID])

    del max_point_cursor        
    del min_point_cursor        

    arcpy.management.SplitLineAtPoint(InputProfiles, max_points, temp_workspace + "\\split_profiles", "1 Meters")
    fieldmappings = arcpy.FieldMappings()
    fieldmappings.addTable(InputProfiles)
    arcpy.SpatialJoin_analysis(temp_workspace + "\\split_profiles", min_points, OutputProfileMetrics, "JOIN_ONE_TO_ONE", "KEEP_COMMON", fieldmappings, "INTERSECT", "1 Meters", "#")
    #arcpy.DeleteField_management(OutputProfileMetrics,['Join_Count','TARGET_FID'])
else:
    arcpy.CopyFeatures_management(InputProfiles, OutputProfileMetrics) 

#a "list" where the name of fields from the attributed table are copied in
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
    

##Axis variables
new_fields = ("profclos", "PI", "profasp","profgrad", "sinuosity") ## All float variables 4
for field in new_fields:
    if field in Fieldlist:
        pass
    else:
        arcpy.AddField_management(OutputProfileMetrics, field, "DOUBLE",10, 2)

if "profamp" in Fieldlist:  ## count = 1
    pass
else:
    #add fieds to attribute tables to be populated with calculated values
    arcpy.AddField_management(OutputProfileMetrics, "profamp", "LONG", 10)


    
##axis curve-fit variables 
new_fields = ("exp_a","exp_b","exp_r2","pow_a", "pow_b", "pow_r2","kcurve_c","kcurve_r2","SL", "SL_r2") ##float variables with high digits
for field in new_fields:
    if field in Fieldlist:
        pass
    else:
        arcpy.AddField_management(OutputProfileMetrics, field, "DOUBLE",10, 4)

##Check the direction and flip the length from low to high elevations
arcpy.AddMessage("Check profile direction and flip it from low to high elevations if necessary...")
Check_If_Flip_Line_Direction(OutputProfileMetrics, InputDEM)

arcpy.AddMessage("Derive profile metrics...")
arcpy.InterpolateShape_3d(InputDEM, OutputProfileMetrics, temp_workspace + "\\profile3D") 

FID_list = []
PI_list = []
HLAsp_list = []
P_clos_list = []
Amplitude_list = []
profgrad_list = []
sinuosity_list = []

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
        if max_Z > min_Z:
            PI = (mean_Z - min_Z) / (max_Z - min_Z)+ 0.001 
        else:
            PI = -999
        #arcpy.AddMessage("High-length HI: " + str(HI))
        PI_list.append(PI)
        Amplitude_list.append(max_Z - min_Z)
        gradient = 180.0/math.pi * math.atan((max_Z - min_Z)/max(LengthfromStart))

        profgrad_list.append(gradient)

        ##derive the sinuosity
        start_end_length = Dist(PointX[0],PointY[0],PointX[-1],PointY[-1])
        #arcpy.AddMessage(start_end_length)
        if start_end_length > 0:
            sinuosity = lineLength / (start_end_length) ##to avoid the division by zero
        else:
            sinuosity = -999
        sinuosity_list.append(sinuosity)

        ##Calculate the HL-Aspect
        #dz  = PointZ[-1] - PointZ[0]
        dx  = PointX[0] - PointX[-1]
        dy  = PointY[0] - PointY[-1]
        #arcpy.AddMessage(str(dx))
        #arcpy.AddMessage(str(dy))

        aspect = 180.0/math.pi * math.atan2(dy, dx)
        #arcpy.AddMessage("Aspect is: " + str(aspect))
        if aspect < 90:
            adj_aspect = 90.0 - aspect
        else:
            adj_aspect = 360 + 90.0 - aspect
        HLAsp_list.append(adj_aspect)

        ##Derive the exponential model fit for the longtitude profile 03/15/2023
        pointH = [y - min_Z for y in PointZ]

        HArr = np.array(pointH)
        LenArr = np.array(LengthfromStart)

        validHArr = HArr[HArr > 0]
        validLenArr = LenArr[HArr > 0]
        
        try:
            polyfit_results = polyfit(validLenArr, np.log(validHArr), 1)
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
        try:
            polyfit_results = polyfit(np.log(validLenArr), np.log(validHArr), 1)
            b = polyfit_results['polynomial'][0]
            a = np.exp(polyfit_results['polynomial'][1])
            R2 = polyfit_results['determination']
        except:
            #arcpy.AddMessage("There is an error!")
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
        #arcpy.AddMessage("the profile closure is: " + str(p_close))
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
        pointH = [y - min_Z for y in PointZ] ##Already reversed in the previous step
        #arcpy.AddMessage(pointH)
        #pointH.reverse()
        #arcpy.AddMessage(pointH)

        #LengthfromStart.reverse() ##Already reversed in the previous step

        pointHArr = np.array(pointH)
        ReverseLengthArr = np.array([(max_len - y) for y in LengthfromStart])

        validHArr = pointHArr[pointHArr > 0]
        #arcpy.AddMessage(validHArr)
        validLenArr = ReverseLengthArr[ReverseLengthArr > 0]
        #arcpy.AddMessage(validLenArr)
        #plt.plot( validLenArr, validHArr)
        #plt.show()
        #plt.plot( np.log(validLenArr), validHArr)
        #plt.show()

        try:
            polyfit_results = polyfit(np.log(validLenArr), validHArr, 1)
            sl = polyfit_results['polynomial'][0]
            #c = np.exp(polyfit_results['polynomial'][1])
            R2 = polyfit_results['determination']
        except:
            #arcpy.AddMessage("There is an error!")
            sl = -999
            #a = -999
            R2 = -999

        SL_list.append(a)
        #SL_c_list.append(b)
        SL_r2_list.append(R2)
        
        i += 1

del row, cursor

fields = ("ProfileID", "profclos", "PI", "profasp", "profamp", "profgrad", "exp_a","exp_b","exp_r2","pow_a", "pow_b", "pow_r2","kcurve_c","kcurve_r2","SL", "SL_r2", "sinuosity")

with arcpy.da.UpdateCursor(OutputProfileMetrics, fields) as cursor:
    for row in cursor:
        try:
            fid = FID_list.index(row[0])
            row[1] = P_clos_list[fid]
            row[2] = PI_list[fid]
            row[3] = HLAsp_list[fid]
            row[4] = Amplitude_list[fid]
            row[5] = profgrad_list[fid]

            row[6] = exp_a_list[fid]
            row[7] = exp_b_list[fid]
            row[8] = exp_r2_list[fid]

            row[9] = pow_a_list[fid]
            row[10] = pow_b_list[fid]
            row[11] = pow_r2_list[fid]

            row[12] = kcurve_c_list[fid]
            row[13] = kcurve_r2_list[fid]
            row[14] = SL_list[fid]
            row[15] = SL_r2_list[fid]
            row[16] = sinuosity_list[fid]

            #update cursor
            cursor.updateRow(row)
        except:
            arcpy.AddMessage("There is an error in the calculation. Move to the next one")
            pass
del row, cursor

arcpy.Delete_management(temp_workspace) ### Empty the in_memory












