import arcpy
import heapq
import tempfile
import math
import sys
import os

# Definitions ---------------------------------------------------------------------------------------------------------
# Priority Queue class
class PQueue:
    def __init__(self):
        self._queue = []
        self._index = 0

    def push(self, item, priority):
        if item not in self._queue:
            heapq.heappush(self._queue, (-priority, self._index, item))
            self._index += 1

    def pop(self):
        self._index -= 1
        return heapq.heappop(self._queue)[-1]

    def length(self):
        return self._index

    def clear(self):
        self._queue = []
        self._index = 0

# Checks outwards from potential summit to make sure it is on the correct spot
def summitBuffer(x,y, buff,raster):

    # Start in top right corner
    x = x + CellSizeX * buff
    y = y + CellSizeY * buff

    # Move down
    for buffer in range(1, (2 * buff) + 1):
        point = str(x) + ' ' + str(y - (buffer * CellSizeY))
        value = arcpy.GetCellValue_management(raster, point).getOutput(0)
        #arcpy.AddMessage(point)
        #arcpy.AddMessage(value)
        if value == 'NoData':
            pass
        else:
            value = float(value)
            neighbors.push(point, value)


    y = y - (CellSizeY * buff * 2)
    # Move right
    for buffer in range(1, (2 * buff) + 1):
        point = str(x - (buffer * CellSizeX)) + ' ' + str(y)
        value = arcpy.GetCellValue_management(raster, point).getOutput(0)
        #arcpy.AddMessage(point)
        #arcpy.AddMessage(value)
        if value == 'NoData':
            pass
        else:
            value = float(value)
            neighbors.push(point, value)


    x = x - (CellSizeX * buff * 2)
    # Move up
    for buffer in range(1, (2 * buff) + 1):
        point = str(x) + ' ' + str(y + (buffer * CellSizeY))
        value = arcpy.GetCellValue_management(raster, point).getOutput(0)
        #arcpy.AddMessage(point)
        #arcpy.AddMessage(value)
        if value == 'NoData':
            pass
        else:
            value = float(value)
            neighbors.push(point, value)


    y = y + (CellSizeY * buff * 2)
    # Move left
    for buffer in range(1, (2 * buff) + 1):
        point = str(x + (buffer * CellSizeX)) + ' ' + str(y)
        value = arcpy.GetCellValue_management(raster, point).getOutput(0)
        #arcpy.AddMessage(point)
        #arcpy.AddMessage(value)
        if value == 'NoData':
            pass
        else:
            value = float(value)
            neighbors.push(point, value)

# Finds the 8 values along with their coords around the xy values sent
def adjValues(x,y,raster):

    # Print inital point to be tested
    point = str(x) + ' ' + str(y)
    print("testing: " + arcpy.GetCellValue_management(raster, point).getOutput(0))

    # Top center
    point = str(x) + ' ' + str(y + CellSizeY)
    value = arcpy.GetCellValue_management(raster, point).getOutput(0)
    print(value)
    if value == 'NoData':
        pass
    else:
        value = float(value)
        neighbors.push(point, value)

    # Top right
    point = str(x + CellSizeX) + ' ' + str(y + CellSizeY)
    value = arcpy.GetCellValue_management(raster, point).getOutput(0)
    #print(value)
    if value == 'NoData':
        pass
    else:
        value = float(value)
        neighbors.push(point, value)

    # Center right
    point = str(x + CellSizeX) + ' ' + str(y)
    value = arcpy.GetCellValue_management(raster, point).getOutput(0)
    #print(value)
    if value == 'NoData':
        pass
    else:
        value = float(value)
        neighbors.push(point, value)

    # Bottom right
    point = str(x + CellSizeX) + ' ' + str(y - CellSizeY)
    value = arcpy.GetCellValue_management(raster, point).getOutput(0)
    #print(value)
    if value == 'NoData':
        pass
    else:
        value = float(value)
        neighbors.push(point, value)

    # Bottom center
    point = str(x) + ' ' + str(y - CellSizeY)
    value = arcpy.GetCellValue_management(raster, point).getOutput(0)
    #print(value)
    if value == 'NoData':
        pass
    else:
        value = float(value)
        neighbors.push(point, value)

    # Bottom left
    point = str(x - CellSizeX) + ' ' + str(y - CellSizeY)
    value = arcpy.GetCellValue_management(raster, point).getOutput(0)
    #print(value)
    if value == 'NoData':
        pass
    else:
        value = float(value)
        neighbors.push(point, value)

    # Center left
    point = str(x - CellSizeX) + ' ' + str(y)
    value = arcpy.GetCellValue_management(raster, point).getOutput(0)
    #print(value)
    if value == 'NoData':
        pass
    else:
        value = float(value)
        neighbors.push(point, value)

    # Top left
    point = str(x - CellSizeX) + ' ' + str(y + CellSizeY)
    value = arcpy.GetCellValue_management(raster, point).getOutput(0)
    #print(value)
    if value == 'NoData':
        pass
    else:
        value = float(value)
        neighbors.push(point, value)

# Helper function to take a string of 2 points and return one of them as float
def seperatePoint(point, index):
    splitPoints = point.split()
    return float(splitPoints[index])

# Variable to keep track of how many pixels a summit moved
global moved
moved = 0

# Checks the neighbors of the point passed to it
def checkNeighbors(testValue, point, raster):
    global summitFound
    while neighbors.length() != 0 and not summitFound:

        poppedPoint = neighbors.pop()
        poppedValue = float(arcpy.GetCellValue_management(raster, poppedPoint).getOutput(0))

        # Value at top of the queue is greater than the value we were testing
        #print("Popped:" + str(poppedValue))
        #print('Test:' + str(testValue))
        if poppedValue > testValue and not summitFound:

            global moved
            moved += 1
            adjValues(seperatePoint(poppedPoint, 0), seperatePoint(poppedPoint, 1), raster)
            checkNeighbors(poppedValue, poppedPoint, raster)

        # Value at top of the queue is not greater than test value should indicate a summit
        if poppedValue <= testValue and not summitFound:

            # Check the summit against the buffer
            for x in range(1, buffer):
                summitBuffer(seperatePoint(point, 0), seperatePoint(point, 1), int(x), raster)

            # Orginal summit holds against buffer
            poppedPoint = neighbors.pop()
            if float(arcpy.GetCellValue_management(raster, poppedPoint).getOutput(0)) <= testValue:

                # The summit was already on the correct pixel
                if moved == 0:

                    print("Summit Found at: " + str(point) + " With a Value of: " + str(testValue))

                    # Find the Distance between the old and new summits
                    distance = 0

                    # Store the data in a dictionary
                    pointsToAdd.append({'NAME': pointName, 'FID': pointFID, 'STATE': pointState, 'COUNTY': pointCounty,
                                        'FEATURE': feature, 'POINT': (seperatePoint(point, 0), seperatePoint(point, 1)),
                                        'DISTANCE': distance, 'ELEV':testValue})
                    # Let the program know a summit has been found, move onto the next point
                    summitFound = True


                # The summit was not on the correct pixel at the start
                else:

                    print("Summit Found at: " + str(point) + " With a Value of: " + str(testValue))

                    # Find the Center of the Pixel
                    colValue = int(abs((xmin - seperatePoint(point, 0)) / CellSizeX) + 1) - 1
                    rowValue = int(abs((ymax - seperatePoint(point, 1)) / CellSizeY) + 1)

                    #print('Row:' + str(rowValue))
                    #print('Col:' + str(colValue))

                    centerX = (xmin + (colValue * CellSizeX)) + (0.5 * CellSizeX)
                    centerY = (ymax - (rowValue * CellSizeY)) + (0.5 * CellSizeY)

                    #print(centerX)
                    #print(centerY)

                    # Find the Distance between the old and new summits
                    a = abs(GNISx - centerX)
                    b = abs(GNISy - centerY)
                    distance = math.sqrt((pow(a, 2)) + (pow(b, 2)))

                    # Store the data in a dictionary
                    pointsToAdd.append({'NAME': pointName, 'FID': pointFID, 'STATE': pointState, 'COUNTY': pointCounty,
                                        'FEATURE': feature, 'POINT': (centerX, centerY), 'DISTANCE': distance, 'ELEV':testValue})
                    # Let the program know a summit has been found, move onto the next point
                    summitFound = True

            # Summit fails against buffer
            else:
                adjValues(seperatePoint(poppedPoint, 0), seperatePoint(poppedPoint, 1), raster)
                checkNeighbors(float(arcpy.GetCellValue_management(raster, poppedPoint).getOutput(0)), poppedPoint, raster)



# Variables -----------------------------------------------------------------------------------------------------------

# Location of the raster image - STR VARIABLE
srcLocation10m = arcpy.GetParameterAsText(0)

# Locations of the shape file to be tested - STR VARIABLE
infc = arcpy.GetParameterAsText(2)

# Buffer for end testing - INT VARIABLE
buffer = int(arcpy.GetParameter(3))

# Directory of the shape file to be created with the results - STR VARIABLE
newFCDir = arcpy.GetParameterAsText(4)

# Name of the shape file to be created with the results - STR VARIABLE
newFCName = arcpy.GetParameterAsText(5)

# Variables (1m)-------------------------------------------------------------------------------------------------------

# Location of the raster image - STR VARIABLE
srcLocation1m = arcpy.GetParameterAsText(1)

# Locations of the shape file to be tested - STR VARIABLE
infc10m = newFCDir + '\\' + newFCName + '.shp'

# Directory of the shape file to be created with the results - STR VARIABLE
newFCDir1m = arcpy.GetParameterAsText(6)

# Name of the shape file to be created with the results - STR VARIABLE
newFCName1m = arcpy.GetParameterAsText(7)

# Temp Workspace location - STR VARIABLE
arcpy.env.workspace = arcpy.GetParameterAsText(8)

# Create a new temporary directory for temp files to be stored in
tempdir = tempfile.mkdtemp()

# Size of the cells in the raster image (10m raster image)
CellSizeX = float(arcpy.GetRasterProperties_management(srcLocation10m, 'CELLSIZEX').getOutput(0))
CellSizeY = float(arcpy.GetRasterProperties_management(srcLocation10m, 'CELLSIZEY').getOutput(0))

# Temp file used with "infc" to get a copy with the same spatial reference as the raster
spatialRef_FC =  tempdir + 'tempFC.shp'

# Temp file created when we convert multipoint fc to point fc
pointFC = tempfile.mkdtemp() + 'tempFC.shp'

# Temp file for clipping the original feature class
clippedFC = tempdir + 'clippedFC.shp'

# Temp directory used to copy the raster
tempSRC = tempdir + 'copiedSRC.shp'

# Priority Queue used to iterate through the neighbors based on their values
neighbors = PQueue()

# Bool used to exit the recursive loop
summitFound = False

# List to be filled with the summit found for each point analyzed
pointsToAdd = []

# Fields to retreive from "infc" and add to the result feature class
#fields = ['gaz_name', 'feature_id', 'state_alph', 'county_nam', 'gaz_featur', 'SHAPE@XY']
fields = ['FEATURE_NA', 'FEATURE_ID', 'STATE_ALPH', 'COUNTY_NAM', 'FEATURE_CL', 'SHAPE@']
#newFields = ['gaz_name', 'feature_id', 'state_alph', 'county_nam', 'gaz_featur', 'SHAPE@XY', 'distance']
newFields = ['FEATURE_NA', 'FEATURE_ID', 'STATE_ALPH', 'COUNTY_NAM', 'FEATURE_CL', 'SHAPE@', 'distance_m', 'elevation']

xmin = float(arcpy.GetRasterProperties_management(srcLocation10m,'LEFT').getOutput(0))
#arcpy.AddMessage(xmin)

ymax = float(arcpy.GetRasterProperties_management(srcLocation10m,'TOP').getOutput(0))
#arcpy.AddMessage(ymax)

# Original X and Y values
global GNISx
global GNISy
GNISx = 0
GNISy = 0

# Code ----------------------------------------------------------------------------------------------------------------
arcpy.env.overwriteOutput = True

# Check that the save directory exists
if os.path.isdir(newFCDir) == False:
    sys.exit('Error: Save Location does not exist!')

# Get spatial analyst license
if arcpy.CheckExtension("Spatial") == "Available":
    arcpy.AddMessage("Checking out Spatial")
    arcpy.CheckOutExtension("Spatial")
else:
    arcpy.AddError("Unable to get spatial analyst extension")
    arcpy.AddMessage(arcpy.GetMessages(0))

# Get the 3D Analysts License
if arcpy.CheckExtension("3D") == "Available":
    arcpy.AddMessage("Checking out 3D")
    arcpy.CheckOutExtension("3D")
else:
    arcpy.AddError("Unable to get 3D analyst extension")
    arcpy.AddMessage(arcpy.GetMessages(0))

# Import sa
from arcpy.sa import *


# Check that the DEM is projected and equidistant (10Meter)
spatial_ref = arcpy.Describe(srcLocation10m).spatialReference
if spatial_ref.factoryCode != 102005 and spatial_ref.factoryCode != 102010:
    sys.exit('Error: DEM must be an equidistant projection! Please use North '
             'America Equidistant Conic or USA Contiguous Equidist conic')

# Check that the DEM is projected and equidistant(1Meter)
spatial_ref = arcpy.Describe(srcLocation1m).spatialReference
if spatial_ref.factoryCode != 102005 and spatial_ref.factoryCode != 102010:
    sys.exit('Error: DEM must be an equidistant projection! Please use North '
             'America Equidistant Conic or USA Contiguous Equidist conic')

# Convert Raster to integer
OutRas = arcpy.sa.Con(Raster(srcLocation10m) > -999, 1, 0)
IntRas = arcpy.Int_3d(OutRas)
arcpy.AddMessage('Raster converted to integer raster')

# Convert the Raster to a polygon feature class
arcpy.RasterToPolygon_conversion(IntRas, tempSRC)
arcpy.AddMessage('Integer raster converted to a polygon feature class')

# Clip the feature class to the extent of the raster
arcpy.Clip_analysis(infc, tempSRC, clippedFC)
arcpy.AddMessage('Feature class clipped to extent of source raster')

# Convert Multipoint feature class to a point feature class
arcpy.MultipartToSinglepart_management(clippedFC, pointFC)
arcpy.AddMessage('Multipoint converted to single point')

# Parse through the points in a feature class, getting rid of mulit-points
# List of feature_id for each row
feature_ids = []
for row in arcpy.da.SearchCursor(pointFC, ["FEATURE_ID", "feature_id"]):
    feature_ids.append(row[0])
multi_points = set([x for n, x in enumerate(feature_ids) if x in feature_ids[:n]])
arcpy.AddMessage('Remaining Multipoint features removed')

# Get the spatial reference of the srcLocation10m
spatialRef = arcpy.Describe(srcLocation10m).spatialReference
arcpy.AddMessage('Spatial reference retrieved')

# Create the Feature Class
arcpy.CreateFeatureclass_management(newFCDir, newFCName, "POINT", infc, "DISABLED", "DISABLED", spatialRef)
#arcpy.AddMessage('New feature class created')

# Reproject the shp
# arcpy.Project_management(pointFC, spatialRef_FC, spatialRef)
# print('Feature class copied from input')

# Find the x,y of GNIS Point
for row in arcpy.da.SearchCursor(pointFC, fields, spatial_reference=spatial_ref):
    pointName = row[0]
    pointFID = row[1]
    pointState = row[2]
    pointCounty = row[3]
    pointShape = row[5]
    x, y = pointShape.firstPoint.X, pointShape.firstPoint.Y
    feature = row[4]
    # If a multipoint, store it, but do not change it
    point = str(x) + ' ' + str(y)
    if row[1] in multi_points and feature == 'Summit':
        pointsToAdd.append({'NAME': pointName, 'FID': pointFID, 'STATE': pointState, 'COUNTY': pointCounty,
                            'FEATURE': feature, 'POINT': pointShape, 'DISTANCE': 0, 'ELEV':float(arcpy.GetCellValue_management(srcLocation10m, point).getOutput(0))})
        continue
    # If a single summit, process through the summit analysis
    if feature == 'Summit':
        arcpy.AddMessage(pointName)
        #arcpy.AddMessage(pointFID)
        #arcpy.AddMessage(feature)
        # Set global bool to false indicating a summit has not been found
        summitFound = False

        # Clear the queue
        neighbors.clear()

        # Set moved to 0 indicating a new summit
        moved = 0

        # Store the original x and y values
        GNISx = x
        GNISy = y

        # Create the point using xy
        point = str(x) + ' ' + str(y)

        # Get value under point
        testValue = float(arcpy.GetCellValue_management(srcLocation10m, point).getOutput(0))

        # Get its 8 neighbors values
        adjValues(seperatePoint(point, 0), seperatePoint(point, 1), srcLocation10m)

        # Test for a summit
        checkNeighbors(testValue, point, srcLocation10m)

    else:
        pass

#arcpy.AddMessage('Writing data to feature class....')
# Add data to the new feature class
FCFile = newFCDir + '\\' + newFCName
if not FCFile.endswith('.shp'):
    FCFile = FCFile + '.shp'
arcpy.AddField_management(FCFile, 'distance_m', 'FLOAT')
arcpy.AddField_management(FCFile, 'elevation', 'FLOAT')
# cursor = arcpy.da.InsertCursor(FCFile, newFields)
for entry in pointsToAdd:
    with arcpy.da.InsertCursor(FCFile, newFields) as cursor:
        row = [entry['NAME'], entry['FID'], entry['STATE'], entry['COUNTY'], entry['FEATURE'], entry['POINT'],
               entry['DISTANCE'], entry['ELEV']]
        cursor.insertRow(row)
        #arcpy.AddMessage(row)
#arcpy.AddMessage('Done!')

pointsToAdd = []
neighbors.clear()




# Code for 1m Raster----------------------------------------------------------------------------------------------------


CellSizeX = float(arcpy.GetRasterProperties_management(srcLocation1m, 'CELLSIZEX').getOutput(0))
CellSizeY = float(arcpy.GetRasterProperties_management(srcLocation1m, 'CELLSIZEY').getOutput(0))

for row in arcpy.da.SearchCursor(infc10m, fields, spatial_reference=spatial_ref):
    pointName = row[0]
    pointFID = row[1]
    pointState = row[2]
    pointCounty = row[3]
    pointShape = row[5]
    x, y = pointShape.firstPoint.X, pointShape.firstPoint.Y
    feature = row[4]
    # If a multipoint, store it, but do not change it
    point = str(x) + ' ' + str(y)
    if row[1] in multi_points and feature == 'Summit':
        pointsToAdd.append({'NAME': pointName, 'FID': pointFID, 'STATE': pointState, 'COUNTY': pointCounty,
                            'FEATURE': feature, 'POINT': pointShape, 'DISTANCE': 0, 'ELEV':float(arcpy.GetCellValue_management(srcLocation1m, point).getOutput(0))})
        continue
    # If a single summit, process through the summit analysis
    if feature == 'Summit':
        arcpy.AddMessage(pointName + ' 1Meter Data')
        #arcpy.AddMessage(pointFID)
        #arcpy.AddMessage(feature)
        # Set global bool to false indicating a summit has not been found
        summitFound = False

        # Clear the queue
        neighbors.clear()

        # Set moved to 0 indicating a new summit
        moved = 0

        # Store the original x and y values
        GNISx = x
        GNISy = y

        # Create the point using xy
        point = str(x) + ' ' + str(y)

        # Get value under point
        testValue = float(arcpy.GetCellValue_management(srcLocation1m, point).getOutput(0))

        # Get its 8 neighbors values
        adjValues(seperatePoint(point, 0), seperatePoint(point, 1), srcLocation1m)

        # Test for a summit
        checkNeighbors(testValue, point, srcLocation1m)

    else:
        pass


# Writing result to file -----------------------------------------------------------------------------------------------
# Create the Feature Class
arcpy.CreateFeatureclass_management(newFCDir1m, newFCName1m, "POINT", infc, "DISABLED", "DISABLED", spatialRef)
arcpy.AddMessage('New feature class created')

arcpy.AddMessage('Writing data to feature class....')
# Add data to the new feature class
FCFile = newFCDir1m + '\\' + newFCName1m
if not FCFile.endswith('.shp'):
    FCFile = FCFile + '.shp'
arcpy.AddField_management(FCFile, 'distance_m', 'FLOAT')
arcpy.AddField_management(FCFile, 'elevation', 'FLOAT')
# cursor = arcpy.da.InsertCursor(FCFile, newFields)
for entry in pointsToAdd:
    with arcpy.da.InsertCursor(FCFile, newFields) as cursor:
        row = [entry['NAME'], entry['FID'], entry['STATE'], entry['COUNTY'], entry['FEATURE'], entry['POINT'],
               entry['DISTANCE'], entry['ELEV']]
        cursor.insertRow(row)
    arcpy.AddMessage(row)
arcpy.AddMessage('Done!')














