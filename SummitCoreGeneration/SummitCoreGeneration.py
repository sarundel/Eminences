# Import sa
from arcpy.sa import *
import heapq
import math
import sys
import os.path
import urllib
import os
import requests
import shutil
import arcpy
import tempfile
import urllib.request
import zipfile
import csv
import glob


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

# Get all neighbors when searching for a possible KeyCol
def GetNeighbors(row, column, queue):
    global MinElevCell
    visited.append((row, column))

    try:
        tl = (row - 1, column + 1)
        value = rasterArray[tl]
        if tl not in visited:
            visited.append(tl)
            if value is None:
                arcpy.AddMessage('NoData')
                pass
            else:
                queue.push(tl, value)
    except Exception as err:
        arcpy.AddMessage(err)

    try:
        tc = (row, column + 1)
        value = rasterArray[tc]
        if tc not in visited:
            visited.append(tc)
            if value is None:
                arcpy.AddMessage('NoData')
                pass
            else:
                queue.push(tc, value)
    except Exception as err:
        arcpy.AddMessage(err)

    try:
        tr = (row + 1, column + 1)
        value = rasterArray[tr]
        if tr not in visited:
            visited.append(tr)
            if value is None:
                arcpy.AddMessage('NoData')
                pass
            else:
                queue.push(tr, value)
    except Exception as err:
        arcpy.AddMessage(err)

    try:
        cr = (row + 1, column)
        value = rasterArray[cr]
        if cr not in visited:
            visited.append(cr)
            if value is None:
                arcpy.AddMessage('NoData')
                pass
            else:
                queue.push(cr, value)
    except Exception as err:
        arcpy.AddMessage(err)

    try:
        br = (row + 1, column - 1)
        value = rasterArray[br]
        if br not in visited:
            visited.append(br)
            if value is None:
                arcpy.AddMessage('NoData')
                pass
            else:
                queue.push(br, value)
    except Exception as err:
        arcpy.AddMessage(err)

    try:
        bc = (row, column - 1)
        value = rasterArray[bc]
        if bc not in visited:
            visited.append(bc)
            if value is None:
                arcpy.AddMessage('NoData')
                pass
            else:
                queue.push(bc, value)
    except Exception as err:
        arcpy.AddMessage(err)

    try:
        bl = (row - 1, column - 1)
        value = rasterArray[bl]
        if bl not in visited:
            visited.append(bl)
            if value is None:
                arcpy.AddMessage('NoData')
                pass
            else:
                queue.push(bl, value)
    except Exception as err:
        arcpy.AddMessage(err)

    try:
        cl = (row - 1, column)
        value = rasterArray[cl]
        if cl not in visited:
            visited.append(cl)
            if value is None:
                arcpy.AddMessage('NoData')
                pass
            else:
                queue.push(cl, value)

    except Exception as err:
        arcpy.AddMessage(err)


def PeakRelDrop(poppedcell, cell_value):
    global height_matched
    global previous
    global MinElevCell

    # If current value is greater than previous value and is greater than original summits value we can stop
    if cell_value > SummitElevValue:
        height_matched = True
    if cell_value < rasterArray[MinElevCell]:
        MinElevCell = poppedcell


# Helper function to take a string of 2 points and return one of them as float
def seperatePoint(point, index):
    splitPoints = point.split()
    return float(splitPoints[index])

# Checks outwards from potential summit to make sure it is on the correct spot
def summitBuffer(x,y, buff):

    # Start in top right corner
    x = x + CellSizeX * buff
    y = y + CellSizeY * buff

    # Move down
    for buffer in range(1, (2 * buff) + 1):
        point = str(x) + ' ' + str(y - (buffer * CellSizeY))
        value = arcpy.GetCellValue_management(srcLocation, point).getOutput(0)
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
        value = arcpy.GetCellValue_management(srcLocation, point).getOutput(0)
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
        value = arcpy.GetCellValue_management(srcLocation, point).getOutput(0)
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
        value = arcpy.GetCellValue_management(srcLocation, point).getOutput(0)
        #arcpy.AddMessage(point)
        #arcpy.AddMessage(value)
        if value == 'NoData':
            pass
        else:
            value = float(value)
            neighbors.push(point, value)

# Finds the 8 values along with their coords around the xy values sent
def adjValues(x,y):

    # Print inital point to be tested
    point = str(x) + ' ' + str(y)
    #print("testing: " + arcpy.GetCellValue_management(srcLocation, point).getOutput(0))

    # Top center
    point = str(x) + ' ' + str(y + CellSizeY)
    value = arcpy.GetCellValue_management(srcLocation, point).getOutput(0)
    #print(value)
    if value == 'NoData':
        pass
    else:
        value = float(value)
        neighbors.push(point, value)

    # Top right
    point = str(x + CellSizeX) + ' ' + str(y + CellSizeY)
    value = arcpy.GetCellValue_management(srcLocation, point).getOutput(0)
    #print(value)
    if value == 'NoData':
        pass
    else:
        value = float(value)
        neighbors.push(point, value)

    # Center right
    point = str(x + CellSizeX) + ' ' + str(y)
    value = arcpy.GetCellValue_management(srcLocation, point).getOutput(0)
    #print(value)
    if value == 'NoData':
        pass
    else:
        value = float(value)
        neighbors.push(point, value)

    # Bottom right
    point = str(x + CellSizeX) + ' ' + str(y - CellSizeY)
    value = arcpy.GetCellValue_management(srcLocation, point).getOutput(0)
    #print(value)
    if value == 'NoData':
        pass
    else:
        value = float(value)
        neighbors.push(point, value)

    # Bottom center
    point = str(x) + ' ' + str(y - CellSizeY)
    value = arcpy.GetCellValue_management(srcLocation, point).getOutput(0)
    #print(value)
    if value == 'NoData':
        pass
    else:
        value = float(value)
        neighbors.push(point, value)

    # Bottom left
    point = str(x - CellSizeX) + ' ' + str(y - CellSizeY)
    value = arcpy.GetCellValue_management(srcLocation, point).getOutput(0)
    #print(value)
    if value == 'NoData':
        pass
    else:
        value = float(value)
        neighbors.push(point, value)

    # Center left
    point = str(x - CellSizeX) + ' ' + str(y)
    value = arcpy.GetCellValue_management(srcLocation, point).getOutput(0)
    #print(value)
    if value == 'NoData':
        pass
    else:
        value = float(value)
        neighbors.push(point, value)

    # Top left
    point = str(x - CellSizeX) + ' ' + str(y + CellSizeY)
    value = arcpy.GetCellValue_management(srcLocation, point).getOutput(0)
    #print(value)
    if value == 'NoData':
        pass
    else:
        value = float(value)
        neighbors.push(point, value)

# Checks the neighbors of the point passed to it
def checkNeighbors(testValue, point, srcLocation,list = []):
    pointName = list[0]
    pointFID = list[1]
    pointState = list[2]
    pointCounty = list[3]
    feature = list[4]

    # update the cell size for testing/calculating the summit points
    CellSizeX = float(arcpy.GetRasterProperties_management(srcLocation, 'CELLSIZEX').getOutput(0))
    CellSizeY = float(arcpy.GetRasterProperties_management(srcLocation, 'CELLSIZEY').getOutput(0))

    global summitFound
    while neighbors.length() != 0 and not summitFound:

        poppedPoint = neighbors.pop()
        poppedValue = float(arcpy.GetCellValue_management(srcLocation, poppedPoint).getOutput(0))

        # Value at top of the queue is greater than the value we were testing
        #print("Popped:" + str(poppedValue))
        #print('Test:' + str(testValue))
        if poppedValue > testValue and not summitFound:

            global moved
            moved += 1
            adjValues(seperatePoint(poppedPoint, 0), seperatePoint(poppedPoint, 1))
            checkNeighbors(poppedValue, poppedPoint, srcLocation, list)

        # Value at top of the queue is not greater than test value should indicate a summit
        if poppedValue <= testValue and not summitFound:

            # Check the summit against the buffer
            for x in range(1, buffer):
                summitBuffer(seperatePoint(point, 0), seperatePoint(point, 1), int(x))

            # Orginal summit holds against buffer
            poppedPoint = neighbors.pop()
            if float(arcpy.GetCellValue_management(srcLocation, poppedPoint).getOutput(0)) <= testValue:
                # The summit was already on the correct pixel
                if moved == 0:

                    print("Summit Found at: " + str(point) + " With a Value of: " + str(testValue) + ' NAEC')

                    # Find the Distance between the old and new summits
                    distance = 0

                    #add point to overall 10m DICT
                    pointsToAddOverall.append({'NAME': pointName, 'FID': pointFID, 'STATE': pointState, 'COUNTY': pointCounty,
                                        'FEATURE': feature, 'POINT': (seperatePoint(point, 0), seperatePoint(point, 1)),
                                        'DISTANCE': distance, 'ELEV':testValue})

                    # Let the program know a summit has been found, move onto the next point
                    summitFound = True


                # The summit was not on the correct pixel at the start
                else:

                    print("Summit Found at: " + str(point) + " With a Value of: " + str(testValue) + ' NAEC')

                    # Find the Center of the Pixel
                    colValue = int(abs((xmin - seperatePoint(point, 0)) / CellSizeX) + 1) - 1
                    rowValue = int(abs((ymax - seperatePoint(point, 1)) / CellSizeY) + 1)

                    #print('Row:' + str(rowValue))
                    #print('Col:' + str(colValue))

                    centerX = (xmin + (colValue * CellSizeX)) + (0.5 * CellSizeX)
                    centerY = (ymax - (rowValue * CellSizeY)) + (0.5 * CellSizeY)

                    # Find the Distance between the old and new summits
                    a = abs(GNISx - centerX)
                    b = abs(GNISy - centerY)
                    distance = math.sqrt((pow(a, 2)) + (pow(b, 2)))


                    # add point to pverall 10m DICT
                    pointsToAddOverall.append({'NAME': pointName, 'FID': pointFID, 'STATE': pointState, 'COUNTY': pointCounty,
                                               'FEATURE': feature, 'POINT': (centerX, centerY),'DISTANCE': distance, 'ELEV': testValue})

                    # Let the program know a summit has been found, move onto the next point

                    summitFound = True

            # Summit fails against buffer
            else:
                adjValues(seperatePoint(poppedPoint, 0), seperatePoint(poppedPoint, 1))
                checkNeighbors(float(arcpy.GetCellValue_management(srcLocation, poppedPoint).getOutput(0)), poppedPoint, srcLocation, list)

# Download the GNIS points, return the path for the selected summit file
def GNIS_Download(save_dir,datasets):
    arcpy.env.overwriteOutput = True
    # List of shapefiles used for the merge
    shape_files = []
    dirpath = save_dir
    fileList = os.listdir(dirpath)
    for FileName in fileList:
        os.remove(dirpath + '\\' + FileName)
    url = "https://viewer.nationalmap.gov/tnmaccess/api/products?"
    parameters = {"datasets": datasets}

    response = requests.get(url, params=parameters)

    data = response.json()
    down_len = len(data['items'])
    print ('%d items returned from TNM API'%down_len)
    print ('Data will be save to:', save_dir)
    for i in range(0,down_len-1):
        url_save = data['items'][i]['downloadURL']
        title = url_save.split('/')
        title = title[-1]
        title = title.replace(" ", "_")
        if (title != 'AllStates.zip' and title != 'NationalFile.zip' and title != 'FM_Features.zip' and title != 'GU_Features.zip'):
            print ('Downloading data from:', url_save)
            urllib.request.urlretrieve(url_save, save_dir + '\\' + title)
            # Unzip the data
            arcpy.AddMessage("Unzipping File...")
            zip_ref = zipfile.ZipFile(save_dir + '\\' + title)
            zip_ref.extractall(save_dir)
            zip_ref.close()
            os.remove(save_dir +'\\'+title)

            arcpy.AddMessage("Converting text files to csv...")
            # Convert .txt files into .csv
            for txt_file in glob.glob(os.path.join(save_dir, '*.txt')):
                with open(txt_file, "r", encoding='UTF-8') as input_file:
                    in_txt = csv.reader(input_file, delimiter='|')
                    filename = os.path.splitext(os.path.basename(txt_file))[0] + '.csv'

                    with open(os.path.join(save_dir, filename), 'w', encoding='UTF-8') as output_file:
                        out_csv = csv.writer(output_file)
                        out_csv.writerows(in_txt)
                os.remove(txt_file)

            arcpy.AddMessage("Converting csv to shapefiles...")
            # Convert csv to shapefile
            # All the shapefiles should save to the individual state folder
            for csv_file in glob.glob(os.path.join(save_dir, '*.csv')):
                print(csv_file)
                temp_name = os.path.splitext(os.path.basename(csv_file))[0]
                arcpy.MakeXYEventLayer_management(csv_file, 'PRIM_LONG_DEC', 'PRIM_LAT_DEC', temp_name)
                arcpy.FeatureClassToShapefile_conversion(temp_name, save_dir)
                os.remove(csv_file)
                print ('Selecting Only summit points...')
                arcpy.Select_analysis(save_dir+'\\'+temp_name+'.shp',save_dir+'\\'+temp_name+'Summit.shp','"FEATURE_CL" = \'Summit\'')
                shape_files.append(save_dir+'\\'+temp_name + 'Summit.shp')
            print ('')
    print ('Merging files...')
    arcpy.Merge_management(shape_files, save_dir+'\\GNIS_AllSummitPoint.shp')
    return(save_dir+'\\GNIS_AllSummitPoint.shp')

#Create feature and pass a bbox to Download data from tnm api
def Create_feature(list = []):
    arcpy.AddMessage('creating feature class for BBox: ')

    spatial_ref = arcpy.SpatialReference(6318)
    arcpy.CreateFeatureclass_management(workspace, 'BBox_temp', "POINT", infc, "DISABLED", "DISABLED",
                                        spatial_ref)
    # Add data to the new feature class
    FCFile = workspace + '\\BBox_temp'
    FCFile_prj = FCFile + "_Project"
    if not FCFile.endswith('.shp'):
        FCFile = FCFile + '.shp'
    if not FCFile_prj.endswith('.shp'):
        FCFile_prj = FCFile_prj + '.shp'
    deleFields_BBox = ['FEATURE_ID', 'STATE_ALPH', 'COUNTY_NAM', 'FEATURE_CL']
    arcpy.DeleteField_management(FCFile, deleFields_BBox)
    # cursor = arcpy.da.InsertCursor(FCFile, newFields)
    for entry in list:
        with arcpy.da.InsertCursor(FCFile, ['FEATURE_NA', 'SHAPE@']) as cursor:
            # 'NAME': pointName,'X-MIN': xmin_bbox, 'Y-MIN': ymin_bbox, 'X-MAX': xmax_bbox, 'Y-MAX': ymax_bbox
            features_name = entry['NAME']
            row = [entry['NAME'], entry['POINT']]
            print (row)
            cursor.insertRow(row)
        #arcpy.AddMessage(row)
    out_corsys = arcpy.SpatialReference(6318)
    #arcpy.AddMessage('AddXY....'
    #arcpy.AddXY_management(FCFile)
    arcpy.AddGeometryAttributes_management(FCFile,'POINT_X_Y_Z_M','KILOMETERS','#',out_corsys)

    arcpy.AddMessage('Converting to different spatial Reference....')

    arcpy.Project_management(FCFile, FCFile_prj, out_corsys)

    # take the xmin, ymax, xmax, ymin from the feature class
    bbox_create = []
    with arcpy.da.SearchCursor(FCFile_prj, ['POINT_X', 'POINT_Y']) as cursor:
        for row in cursor:
            bbox_create.append(row[0])
            bbox_create.append(row[1])
    # swap the position of ymax and ymin to fit the API
    temp_swap = bbox_create[1]
    bbox_create[1] = bbox_create[3]
    bbox_create[3] = temp_swap

    prodExtents = "1 x 1 degree"
    datasets = "National Elevation Dataset (NED) 1/3 arc-second"
    arcpy.AddMessage('Requesting to get files from TNM....')

    TNM_API_Download(features_name,datasets, bbox_create, prodExtents)

    # we dont need the BBox.shp file for the next point. Delete it to save space
    arcpy.DeleteFeatures_management(FCFile)
    spat_ref_NAEC = arcpy.SpatialReference('North America Equidistant Conic')
    if len(merge_file) > 1:
        arcpy.AddMessage("Merge Raster...")
        arcpy.MosaicToNewRaster_management(merge_file, workspace, '\\MergedIMG_NAEC.tif',
                                           coordinate_system_for_the_raster=spat_ref_NAEC,
                                           pixel_type='32_BIT_FLOAT',
                                           number_of_bands=1, mosaic_method='LAST', mosaic_colormap_mode='FIRST')
        return (workspace + '\\MergedIMG_NAEC.tif')
    elif len(merge_file) == 0:
        print ('No Raster IMG found...')
    else:
        arcpy.MosaicToNewRaster_management(merge_file[0], workspace, '\\MergedIMG_NAEC.tif',
                                           coordinate_system_for_the_raster=spat_ref_NAEC,
                                           pixel_type='32_BIT_FLOAT',
                                           number_of_bands=1, mosaic_method='LAST', mosaic_colormap_mode='FIRST')
        return (workspace + '\\MergedIMG_NAEC.tif')

#Download raster from TNM and return the files path
def TNM_API_Download(Feature_name, datasets, bbox, prodExtents):
    arcpy.env.overwriteOutput = True
    url = "https://viewer.nationalmap.gov/tnmaccess/api/products?"
    parameters = {"bbox": bbox, "datasets": datasets, "prodFormats": 'IMG', "prodExtents": prodExtents}

    response = requests.get(url, params=parameters)

    data = response.json()
    down_len = len(data['items'])
    print ('%d items returned from TNM API'%down_len)
    print ('Data will be save to:', workspace)
    if down_len == 0:
        print ('No '+ datasets +' for ' + Feature_name+'...')
        delete_feature.append(Feature_name)
    else:
        for i in range(0,down_len):
            url_save = data['items'][i]['downloadURL']
            print("Downloading info from requested web...")
            title = url_save.split('/')
            title = title[-1]
            title = title.replace(" ", "_")
            file_name = title[:-4] + '.img'
            if len(merge_file) == 0:
                merge_file.append(workspace + '\\' + title[:-4] + '\\' + file_name)
            else:
                if merge_file[0] != workspace + '\\' + title[:-4] + '\\' + file_name:
                    merge_file.insert(0,workspace + '\\' + title[:-4] + '\\' + file_name)

            if os.path.exists(workspace + '\\' + title):
                print (title + ': Already exists')
            else:
                print ('Downloading data from:', url_save)
                urllib.request.urlretrieve(url_save, workspace + '\\' + title)
                # Unzip the data
                arcpy.AddMessage("Unzipping File...")
                zip_ref = zipfile.ZipFile(workspace + '\\' + title)
                title = title[:-4]
                os.makedirs(workspace + '\\' + title)
                workspace_zip = workspace + '\\' + title
                zip_ref.extractall(workspace_zip)
                zip_ref.close()

# Find the summit points
def FindSummit(source_location, input_file, store_dir, store_name):

    arcpy.CreateFeatureclass_management(store_dir, store_name, "POINT", input_file, "DISABLED", "DISABLED", spatial_ref)

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
        point_list = [pointName, pointFID, pointState, pointCounty, feature, pointShape]

        if row[1] in multi_points and feature == 'Summit':
            pointsToAddOverall.append({'NAME': pointName, 'FID': pointFID, 'STATE': pointState, 'COUNTY': pointCounty,
                                'FEATURE': feature, 'POINT': pointShape, 'DISTANCE': 0,
                                'ELEV': float(arcpy.GetCellValue_management(source_location, point).getOutput(0))})
            continue

        # If a single summit, process through the summit analysis
        if feature == 'Summit':
            print ('')
            arcpy.AddMessage(pointName)
            arcpy.AddMessage(pointFID)
            arcpy.AddMessage(feature)
            # Set global bool to false indicating a summit has not been found
            global summitFound
            summitFound = False

            # Clear the queue
            neighbors.clear()

            global moved
            # Set moved to 0 indicating a new summit
            moved = 0

            global GNISx,GNISy
            # Store the original x and y values
            GNISx = x
            GNISy = y

            # Create the point using xy
            point = str(x) + ' ' + str(y)

            # Get value under point
            testValue = float(arcpy.GetCellValue_management(source_location, point).getOutput(0))

            # Get its 8 neighbors values
            adjValues(seperatePoint(point, 0), seperatePoint(point, 1))

            # Test for a summit
            checkNeighbors(testValue, point, source_location, point_list)

        else:
            pass

    arcpy.AddMessage('Writing data to feature class....')
    # Add data to the new feature class
    FCFile = store_dir + '\\' + store_name
    if not FCFile.endswith('.shp'):
        FCFile = FCFile + '.shp'
    arcpy.AddField_management(FCFile, 'distance_m', 'FLOAT')
    arcpy.AddField_management(FCFile, 'elevation', 'FLOAT')
    # cursor = arcpy.da.InsertCursor(FCFile, newFields)
    for entry in pointsToAddOverall:
        with arcpy.da.InsertCursor(FCFile, newFields) as cursor:
            row = [entry['NAME'], entry['FID'], entry['STATE'], entry['COUNTY'], entry['FEATURE'], entry['POINT'],
                    entry['DISTANCE'], entry['ELEV']]
            cursor.insertRow(row)
            arcpy.AddMessage(row)
    arcpy.AddMessage('Done!')

    return (FCFile)

'''-----------------------------
----------------------------------'''
# Variables for the user define functions

# Variable to keep track of how many pixels a summit moved
global moved
moved = 0
# temp[] is use to make the polygon to clip the shapefile
temp = []
temp2 = []
# bbox[] is use to download/get the srcLocation (TIF/IMG)
bbox = []
merge_file = []
# Bool used to exit the recursive loop
summitFound = False
# List to be filled with the summit found for each point analyzed
pointsToAddOverall = []
# List for creating new feature class
pointsToAdd = []
# final feature class
pointsfinals = []
# delete the feature that dont have raster data
delete_feature = []
# Fields to retreive from "infc" and add to the result feature class
#fields = ['gaz_name', 'feature_id', 'state_alph', 'county_nam', 'gaz_featur', 'SHAPE@XY']
fields = ['FEATURE_NA', 'FEATURE_ID', 'STATE_ALPH', 'COUNTY_NAM', 'FEATURE_CL', 'SHAPE@']
#newFields = ['gaz_name', 'feature_id', 'state_alph', 'county_nam', 'gaz_featur', 'SHAPE@XY', 'distance']
newFields = ['FEATURE_NA', 'FEATURE_ID', 'STATE_ALPH', 'COUNTY_NAM', 'FEATURE_CL', 'SHAPE@', 'distance_m', 'elevation']
# Create a new temporary directory for temp files to be stored in
tempdir = tempfile.mkdtemp()
# Temp file used with "infc" to get a copy with the same spatial reference as the raster
spatialRef_FC =  tempdir + 'tempFC.shp'
# Temp file created when we convert multipoint fc to point fc
pointFC = tempfile.mkdtemp() + 'tempFC.shp'
# Temp file for clipping the original feature class
clippedFC = tempdir + 'clippedFC.shp'
# Temp directory used to copy the raster
tempSRC = tempdir + 'copiedSRC.shp'
# Buffer for END testing, should be an INT
buffer = 7
# Original X and Y values
global GNISx
global GNISy
GNISx = 0
GNISy = 0

# A checkbox to ask the user whether they want to download gnis_data
Checkbox_GNIS = arcpy.GetParameterAsText(0)

# Dir that would be saved for the GNIS shapefiles
Save_Dir_GNIS = arcpy.GetParameterAsText(1)

# Bounding box to creat the TIF files and polygon to clip the shapefile(infc) (xin,ymin,xmax,ymax) in NAEC
BBox_Value = arcpy.GetParameterAsText(2)

#Dir that will save the core generation result
Dir = arcpy.GetParameterAsText(3)

name_save = arcpy.GetParameterAsText(4)

# Ask user for the infc, shapfile if they don't want to download the update
infc = arcpy.GetParameterAsText(5)

# Bool to decide if user wants polygons
out_poly = 'true'

# Workspace, all the files that the program generate will be store there
workspace = arcpy.GetParameterAsText(6)

arcpy.env.workspace = workspace
# Priority Queue used to iterate through the neighbors based on their values
neighbors = PQueue()


'''--------------------------------Main------------------------------------'''
arcpy.env.overwriteOutput = True


# download the gnis data
if str(Checkbox_GNIS) == 'true':
    # download the gnis_data, infc = input feature class, the shapefile(summit)
    infc = GNIS_Download(Save_Dir_GNIS,'National Geographic Names Information System (GNIS)')
    temp = BBox_Value.split(',')
    for entries in temp:
        temp2.append(float(entries))
    # create an arrat of point to make the polygon
    array = arcpy.Array([arcpy.Point(temp[0], temp[3]), arcpy.Point(temp[2], temp[1]), arcpy.Point(temp[0], temp[1]),
                         arcpy.Point(temp[2], temp[3])])
    # specific format for the create_feature function
    bbox.append({'NAME': 'GNIS_BBOX', 'POINT': (temp2[0], temp2[3])})
    bbox.append({'NAME': 'GNIS_BBOX', 'POINT': (temp2[2], temp2[1])})

    # make a poly to clip gnis points
    polygon = arcpy.Polygon(array, arcpy.SpatialReference(6318))
    # clip the GNIS Shape file according to the bbox
    arcpy.Clip_analysis(infc,polygon,infc[:-4]+'_StudyArea.shp')
    infc = infc[:-4]+'_StudyArea.shp'
else:
    arcpy.AddMessage("Selecting only summit points")
    arcpy.Select_analysis(infc, infc[:-4]+'_Summit.shp',
                          '"FEATURE_CL" = \'Summit\'')
    infc = infc[:-4]+'_Summit.shp'

    arcpy.MinimumBoundingGeometry_management(infc, workspace+'\\poly.shp', geometry_type='RECTANGLE_BY_WIDTH')
    out_cor = arcpy.SpatialReference(6318)
    arcpy.AddGeometryAttributes_management(workspace+'\\poly.shp',Geometry_Properties='EXTENT',Coordinate_System=out_cor)
    with arcpy.da.SearchCursor(workspace+'\\poly.shp', ['EXT_MIN_X', 'EXT_MIN_Y','EXT_MAX_X','EXT_MAX_Y']) as cursor:
        for row in cursor:
            bbox.append({'NAME': 'GNIS_BBOX', 'POINT': (float(row[0]),float(row[3]))})
            bbox.append({'NAME': 'GNIS_BBOX', 'POINT': (float(row[2]), float(row[1]))})


# Create_feature will return the path for the clipped TIF/IMG
srcLocation = Create_feature(bbox)
merge_file = []
xmin = float(arcpy.GetRasterProperties_management(srcLocation,'LEFT').getOutput(0))
ymax = float(arcpy.GetRasterProperties_management(srcLocation,'TOP').getOutput(0))
# Size of the cells in the raster image
CellSizeX = float(arcpy.GetRasterProperties_management(srcLocation, 'CELLSIZEX').getOutput(0))
CellSizeY = float(arcpy.GetRasterProperties_management(srcLocation, 'CELLSIZEY').getOutput(0))
# get the spatial_ref from the source location, should be North America Equidistant Coic
spatial_ref = arcpy.Describe(srcLocation).spatialReference


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

# Check that the DEM is projected and equidistant (10Meter)
spatial_ref = arcpy.Describe(srcLocation).spatialReference
if spatial_ref.factoryCode != 102005 and spatial_ref.factoryCode != 102010:
    sys.exit('Error: DEM must be an equidistant projection! Please use North '
             'America Equidistant Conic or USA Contiguous Equidist conic')

# Convert Raster to integer
OutRas = arcpy.sa.Con(Raster(srcLocation) > -999, 1, 0)
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


# Find the summit points for 1/3(10m) Raster, return the path to the shapefile(infc)
TenMeterSummit = FindSummit(srcLocation, infc, workspace,'TenMeterSummit')

DEM = srcLocation
infc = TenMeterSummit

# Temp file for clipping the original feature class
clippedFC_Poly = tempdir + 'clippedFC_Poly.shp'

# Temp directory used to copy the raster
tempDEM_Poly = tempdir + 'copiedSRC_Poly.shp'

# Largest value in the grid, indicating the highest cell
MAX = arcpy.GetRasterProperties_management(DEM, "MAXIMUM").getOutput(0)

# Fields for the searchcursor (two different ones, for different shape files)
fields_Poly = ['SHAPE@', 'FEATURE_NA', 'FEATURE_ID', 'STATE_ALPH', 'COUNTY_NAM', 'FEATURE_CL']


# Raster Properties for geting row & column of a pixel
xmin = float(arcpy.GetRasterProperties_management(DEM, 'LEFT').getOutput(0))
ymax = float(arcpy.GetRasterProperties_management(DEM, 'TOP').getOutput(0))
ymin = float(arcpy.GetRasterProperties_management(DEM, 'BOTTOM').getOutput(0))

pixelSizeX = float(arcpy.GetRasterProperties_management(DEM, 'CELLSIZEX').getOutput(0))
pixelSizeY = float(arcpy.GetRasterProperties_management(DEM, 'CELLSIZEY').getOutput(0))

lowerLeft = arcpy.Point(xmin, ymin)

spatialRef = arcpy.Describe(infc).spatialReference
# Code ----------------------------------------------------------------------------------------------------------------

arcpy.env.overwriteOutput = True

arcpy.env.scratchWorkspace = workspace

arcpy.env.outputCoordinateSystem = arcpy.Describe(DEM).spatialReference

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

# Convert Raster to integer
OutRas = arcpy.sa.Con(Raster(DEM) > -999, 1, 0)
IntRas = arcpy.Int_3d(OutRas)
arcpy.AddMessage('Raster converted to integer raster')

# Convert the Raster to a polygon feature class
arcpy.RasterToPolygon_conversion(IntRas, tempDEM_Poly)
arcpy.AddMessage('Integer raster converted to a polygon feature class')

# Clip the feature class to the extent of the raster
arcpy.Clip_analysis(infc, tempDEM_Poly, clippedFC_Poly)
arcpy.AddMessage('Feature class clipped to extent of source raster')

# Create our queues
neighbors = PQueue()
summitTest = PQueue()

# Skip multi-point feature classes
feature_ids = []
for row in arcpy.da.SearchCursor(clippedFC_Poly, ["FEATURE_ID", "feature_id"]):
    feature_ids.append(row[0])
multi_points = set([x for n, x in enumerate(feature_ids) if x in feature_ids[:n]])

if not name_save.endswith('.shp'):
    name_save += '.shp'

if out_poly == 'true':
    # Create our feature class
    arcpy.AddMessage('Creating polygon feature class')
    arcpy.CreateFeatureclass_management(Dir, name_save, "POLYGON", '#', "DISABLED", "DISABLED", spatialRef)
    polygonfc = Dir + '\\' + name_save
    arcpy.AddField_management(polygonfc, "Name", "TEXT")
    arcpy.AddField_management(polygonfc, "FeatureID", "TEXT")

# Loop through the feature class
for row in arcpy.da.SearchCursor(clippedFC_Poly, fields_Poly):
    # Get the X and Y values
    pointShape = row[0]
    pointX, pointY = pointShape.firstPoint.X, pointShape.firstPoint.Y
    # Get the name of the point
    pointName = row[1]
    # Get the ID of the point
    pointFID = row[2]
    # Get the state the point is located in
    pointState = row[3]
    # Get the county the point is located in
    pointCounty = row[4]
    # Type of feature the point is
    feature = row[5]

    if row[2] in multi_points:
        continue

    # Convert the Raster to a numpy array
    global rasterArray
    rasterArray = arcpy.RasterToNumPyArray(DEM)

    # Get the row & column values from the X & Y values
    global SummitColValue
    global SummitRowValue
    global SummitElevValue
    SummitColValue = int((abs(xmin - pointX)) / pixelSizeX)
    SummitRowValue = int((abs(ymax - pointY)) / pixelSizeY)
    SummitElevValue = rasterArray[SummitRowValue, SummitColValue]

    # Skip the highest summit in the study area
    max = round(float(MAX), 2)
    if (float(max) - float(SummitElevValue)) < .1:
        arcpy.AddMessage("Skipping Max")
        continue

    # Clean up previous iterations
    neighbors.clear()
    height_matched = False
    visited = []

    # Used to keep track of the lowest elevation cell encountered
    global MinElevCell
    MinElevCell = (SummitRowValue, SummitColValue)
    KeyColLocation = MinElevCell

    arcpy.AddMessage("Testing: " + str(pointName))
    # Get the neighbors of the summit
    GetNeighbors(SummitRowValue, SummitColValue, neighbors)
    # Set summit value to previous
    previous = rasterArray[SummitRowValue, SummitColValue]
    # Pop the next highest cell
    poppedcell = neighbors.pop()
    # Check cell for possible KeyCol
    PeakRelDrop(poppedcell, rasterArray[poppedcell])

    while height_matched == False:
        try:

            # Get neighbors of the popped cell
            GetNeighbors(poppedcell[0], poppedcell[1], neighbors)
            # Pop the next highest cell
            poppedcell = neighbors.pop()
            # Check for KeyCol
            PeakRelDrop(poppedcell, rasterArray[poppedcell])
            # Update Key Col
            KeyColLocation = MinElevCell
        except Exception as err:
            arcpy.AddMessage(err)
            # Clear the data structures
    visited = []

    neighbors.clear()
    # Global list of core cells
    global CoreCellLocations
    CoreCellLocations = []

    poppedcell = (SummitRowValue, SummitColValue)
    arcpy.AddMessage("Key Col Found")

    CoreCellLocations.append(poppedcell)

    arcpy.AddMessage("Finding Core")
    KeyColMatched = False

    if poppedcell == KeyColLocation:
        KeyColMatched = True

    while KeyColMatched == False:
        try:
            GetNeighbors(poppedcell[0], poppedcell[1], neighbors)

            CoreCellLocations.append(poppedcell)

            poppedcell = neighbors.pop()

            if poppedcell == KeyColLocation:
                KeyColMatched = True

        except Exception as err:
            arcpy.AddMessage(err)

    CoreCellLocations.append(KeyColLocation)

    for x in range(0, rasterArray.shape[0]):
        for y in range(0, rasterArray.shape[1]):
            rasterArray[x, y] = None

    arcpy.AddMessage("Setting Cell Values")
    for i in range(0, len(CoreCellLocations)):
        rasterArray[CoreCellLocations[i]] = pointFID

    tempRaster = os.path.join(Dir, pointName + '_' + str(pointFID) + '.tif')
    tempRaster = tempRaster.replace(" ", "_")


    arcpy.AddMessage("Converting to Raster")
    numpyRaster = arcpy.NumPyArrayToRaster(rasterArray, lowerLeft, pixelSizeX, pixelSizeY)

    RasterLocation = os.path.join(Dir, pointName + '_Raster' + '.tif')

    arcpy.AddMessage("Saving Raster")
    numpyRaster.save(RasterLocation)

    if out_poly == 'true':

        arcpy.AddMessage("Converting to Int Raster for polygon generation")
        arcpy.Int_3d(numpyRaster, tempRaster)

        fc = Dir + '\\' + pointName + '.shp'
        fc = fc.replace(" ", "_")

        arcpy.AddMessage("Converting to Polygon")
        arcpy.RasterToPolygon_conversion(tempRaster, fc, 'NO_SIMPLIFY', 'Value', 'MULTIPLE_OUTER_PART')

        arcpy.AddField_management(fc, "Name", "TEXT")
        arcpy.AddField_management(fc, "FeatureID", "TEXT")

        arcpy.DeleteField_management(fc, "gridcode")

        with arcpy.da.UpdateCursor(fc, ['Name', 'FeatureID']) as cursor:
            for row in cursor:
                row = [pointName, pointFID]
                print(row)
                cursor.updateRow(row)

        print(polygonfc)
        print(fc)

        arcpy.AddMessage("Appending to feature class")
        arcpy.Append_management(fc, polygonfc, 'TEST')

        arcpy.Delete_management(fc)
