'''
SummitSpotHeightCalculator:
    Will ask for:
        1. GNIS Download
        2. Bounding Box to get the 1/3 data
'''
#Import sa
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
import xml.etree.ElementTree as ET

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

def create_epsg(filename,file):
    os.chdir(workspace+'\\'+filename)
    file = file[:-4] + '.xml'
    print(file)
    print(os.getcwd())

    tree = ET.parse(file)
    root = tree.getroot()
    for auth in root.iter('MapProjectionDefinition'):
        epsg = ''
        start = auth.text.find('EPSG') + 6
        epsg_code = auth.text[start:]
        if epsg_code.isdigit() == False:
            for dig in epsg_code:
                if dig.isdigit():
                    epsg += dig
    wkt = urllib.request.urlopen("http://spatialreference.org/ref/epsg/{0}/prettywkt/".format(epsg))
    remove_spaces = wkt.read().decode('utf-8').replace(' ', '')
    output = remove_spaces.replace('\n', '')
    prj = open(filename + ".prj", 'w')
    prj.write(output)
    prj.close()


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
def checkNeighbors(testValue, point, srcLocation,list):
    pointName = list[0]
    pointFID = list[1]
    pointState = list[2]
    pointCounty = list[3]
    feature = list[4]
    BboxToAdd = []
    # update the cell size for testing/calculating the summit points
    CellSizeX = float(arcpy.GetRasterProperties_management(srcLocation, 'CELLSIZEX').getOutput(0))
    CellSizeY = float(arcpy.GetRasterProperties_management(srcLocation, 'CELLSIZEY').getOutput(0))
    if CellSizeX < 1.1:
        test_data = "onemeter"
        source = '1m'
    elif CellSizeX > 2:
        test_data = "tenmeter"
        source = '1/3'
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

                    if test_data == "onemeter":
                        xmin_bbox = (seperatePoint(point, 0) - CellSizeX * 7.0)
                        ymin_bbox = (seperatePoint(point, 1) - CellSizeY * 7.0)
                        xmax_bbox = (seperatePoint(point, 0) + CellSizeX * 7.0)
                        ymax_bbox = (seperatePoint(point, 1) + CellSizeY * 7.0)
                    else:
                        xmin_bbox = (seperatePoint(point, 0) - CellSizeX * 4.0)
                        ymin_bbox = (seperatePoint(point, 1) - CellSizeY * 4.0)
                        xmax_bbox = (seperatePoint(point, 0) + CellSizeX * 4.0)
                        ymax_bbox = (seperatePoint(point, 1) + CellSizeY * 4.0)

                    # Find the Distance between the old and new summits
                    distance = 0

                    #add point to overall 10m DICT
                    pointsToAddOverall.append({'NAME': pointName, 'FID': pointFID, 'STATE': pointState, 'COUNTY': pointCounty,
                                        'FEATURE': feature, 'POINT': (seperatePoint(point, 0), seperatePoint(point, 1)),
                                        'DISTANCE': distance, 'ELEV':testValue, 'SOURCE': source})

                    # Store the bbox data in a dictionary
                    BboxToAdd.append({'NAME': pointName, 'POINT': (xmin_bbox, ymax_bbox)})
                    BboxToAdd.append({'NAME': pointName, 'POINT': (xmax_bbox, ymin_bbox)})

                    #Create a feature class for the value and download the raster that we need
                    Create_feature(test_data,BboxToAdd)

                    # Let the program know a summit has been found, move onto the next point
                    summitFound = True
                    BboxToAdd = []


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

                    if test_data == "onemeter":
                        xmin_bbox = (centerX - CellSizeX * 7.0)
                        ymin_bbox = (centerY - CellSizeY * 7.0)
                        xmax_bbox = (centerX + CellSizeX * 7.0)
                        ymax_bbox = (centerY + CellSizeY * 7.0)
                    else:
                        xmin_bbox = (centerX - CellSizeX * 4.0)
                        ymin_bbox = (centerY - CellSizeY * 4.0)
                        xmax_bbox = (centerX + CellSizeX * 4.0)
                        ymax_bbox = (centerY + CellSizeY * 4.0)

                    # Find the Distance between the old and new summits
                    a = abs(GNISx - centerX)
                    b = abs(GNISy - centerY)
                    distance = math.sqrt((pow(a, 2)) + (pow(b, 2)))


                    # add point to pverall 10m DICT
                    pointsToAddOverall.append({'NAME': pointName, 'FID': pointFID, 'STATE': pointState, 'COUNTY': pointCounty,
                                               'FEATURE': feature, 'POINT': (centerX, centerY),'DISTANCE': distance, 'ELEV': testValue,'SOURCE': source})

                    BboxToAdd.append({'NAME': pointName, 'POINT': (xmin_bbox, ymax_bbox)})
                    BboxToAdd.append({'NAME': pointName, 'POINT': (xmax_bbox, ymin_bbox)})

                    #Create a feature class for the value and download the raster that we need
                    Create_feature(test_data,BboxToAdd)

                    # Let the program know a summit has been found, move onto the next point

                    summitFound = True
                    BboxToAdd = []

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
def Create_feature(TestData, list):
    arcpy.AddMessage('creating feature class for BBox: ')
    if TestData != 'GNIS':
        spatial_ref = arcpy.Describe(srcLocation).spatialReference
    else:
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

    arcpy.AddMessage('Requesting to get files from TNM....')

    # take the xmin, ymax, xmax, ymin from the feature class
    bbox_create = []
    with arcpy.da.SearchCursor(FCFile_prj, ['POINT_X', 'POINT_Y']) as cursor:
        for row in cursor:
            bbox_create.append(round(row[0],3))
            bbox_create.append(round(row[1],3))
    # swap the position of ymax and ymin to fit the API
    temp_swap = bbox_create[1]
    bbox_create[1] = bbox_create[3]
    bbox_create[3] = temp_swap
    #dealing with single point
    if bbox_create[0] == bbox_create[2]:
        bbox_create[1] = bbox_create[3]
    elif bbox_create[1] == bbox_create[3]:
        bbox_create[0] = bbox_create[2]

    arcpy.AddMessage(bbox_create)
    # Have to check what kind of data we want
    # 10m feature class -> Dataset: Digital Elevation Model (DEM) 1 meter, prodExtents: 10000 x 10000 meter
    # prodFormats: IMG
    # 1m feature class -> Dataset: Lidar Point Cloud (LPC)
    if TestData == "tenmeter":
        prodExtents = "10000 x 10000 meter"
        datasets = "Digital Elevation Model (DEM) 1 meter"
    elif TestData == "onemeter":
        prodExtents = "Varies"
        datasets = "Lidar Point Cloud (LPC)"
    elif TestData == "GNIS":
        prodExtents = "1 x 1 degree"
        datasets = "National Elevation Dataset (NED) 1/3 arc-second"

    TNM_API_Download(features_name,datasets, bbox_create, prodExtents)

    # we dont need the BBox.shp file for the next point. Delete it to save space
    arcpy.DeleteFeatures_management(FCFile)
    if TestData == 'GNIS':
        spat_ref_NAEC = arcpy.SpatialReference('North America Equidistant Conic')
        if len(merge_file) > 1:
            arcpy.AddMessage("Merging Raster...")
            arcpy.MosaicToNewRaster_management(merge_file, workspace, '\\MergedIMG_NAEC.tif',
                                               coordinate_system_for_the_raster=spat_ref_NAEC, pixel_type='32_BIT_FLOAT',
                                               number_of_bands=1, mosaic_method='LAST', mosaic_colormap_mode='FIRST')
            return (workspace + '\\MergedIMG_NAEC.tif')
        elif len(merge_file) == 0:
            arcpy.AddMessage('No Raster IMG found...')
        else:
            arcpy.MosaicToNewRaster_management(merge_file[0], workspace, '\\MergedIMG_NAEC.tif',
                                               coordinate_system_for_the_raster=spat_ref_NAEC, pixel_type='32_BIT_FLOAT',
                                               number_of_bands=1, mosaic_method='LAST', mosaic_colormap_mode='FIRST')
            return (workspace + '\\MergedIMG_NAEC.tif')

#Download raster from TNM and return the files path
def TNM_API_Download(Feature_name, datasets, bbox, prodExtents):
    arcpy.env.overwriteOutput = True
    las_merge = []
    url = "https://viewer.nationalmap.gov/tnmaccess/api/products?"
    parameters = {"bbox": bbox,"datasets": datasets, "prodExtents": prodExtents}
    if datasets != 'Lidar Point Cloud (LPC)':
        parameters = {"bbox": bbox, "datasets": datasets, "prodFormats": 'IMG', "prodExtents": prodExtents}
        print(bbox, datasets, prodExtents)
    arcpy.AddMessage('Dataset:%s' %datasets)
    arcpy.AddMessage('prodExtents:%s' % prodExtents)
    arcpy.AddMessage('sending request to TNM...')
    response = requests.get(url, params=parameters)
    data = response.json()
    down_len = len(data['items'])
    arcpy.AddMessage('%d items returned from TNM API'%down_len)
    if down_len == 0:
        print ('No '+ datasets +' for ' + Feature_name+'...')
        delete_feature.append(Feature_name)
        print(delete_feature)
    else:
        if datasets != 'Lidar Point Cloud (LPC)':
            for i in range(0,down_len):
                url_save = data['items'][i]['downloadURL']
                title = url_save.split('/')
                title = title[-1]
                title = title.replace(" ", "_")

                if os.path.exists(workspace + '\\' + title):
                    workspace_zip = workspace + '\\' + title[:-4]
                    print(workspace_zip)
                    print (title + ': Already exists')
                    if datasets == 'Lidar Point Cloud (LPC)':
                        file_name = workspace_zip + '\\' + title[:-4] + '.las'
                        if title[0] == 'N':
                            create_epsg(title[:-4], file_name)
                            arcpy.env.workspace = workspace
                        print(file_name)
                        las_merge.append(file_name)
                    if (glob.glob(workspace_zip + '\\*.img')) in merge_file:
                        print(title[:-4] + ' already in list')
                    else:
                        if len(merge_file) == 0:
                            merge_file.append(glob.glob(workspace_zip + '\\*.img'))
                            print('adding file to list')
                        else:
                            merge_file.insert(0, glob.glob(workspace_zip + '\\*.img'))
                            print('insert file to list')

                else:
                    arcpy.AddMessage('Downloading data from %s'%url_save)
                    workspace_zip = workspace + '\\' + title[:-4]
                    urllib.request.urlretrieve(url_save, workspace + '\\' + title)
                    # Unzip the data
                    arcpy.AddMessage("Unzipping File...")
                    if os.path.exists(workspace_zip):
                        arcpy.AddMessage('File already exist...')
                    else:
                        zip_ref = zipfile.ZipFile(workspace + '\\' + title)
                        title = title[:-4]
                        os.makedirs(workspace + '\\' + title)
                        workspace_zip = workspace + '\\' + title
                        print('workspace:',workspace)
                        print('workspace_zip',workspace_zip)
                        zip_ref.extractall(workspace_zip)
                        if datasets == 'Lidar Point Cloud (LPC)':
                            file_name = workspace_zip + '\\' + title + '.las'
                            if title[0] == 'N':
                                create_epsg(title[:-4], file_name)
                                arcpy.env.workspace = workspace
                            las_merge.append(file_name)
                            print(las_merge)
                        if len(merge_file) == 0:
                            merge_file.append(glob.glob(workspace_zip+'\\*.img'))
                        else:
                            merge_file.insert(0, glob.glob(workspace_zip+'\\*.img'))
                        zip_ref.close()
        else:
            url_save = data['items'][0]['downloadURL']
            arcpy.AddMessage("Downloading info from requested web...")
            title = url_save.split('/')
            title = title[-1]
            title = title.replace(" ", "_")

            if os.path.exists(workspace + '\\' + title):
                workspace_zip = workspace + '\\' + title[:-4]
                print(workspace_zip)
                arcpy.AddMessage(title + ': Already exists')
                file_name = workspace_zip + '\\' + title[:-4] + '.las'
                if title[0] == 'N':
                    create_epsg(title[:-4], file_name)
                    arcpy.env.workspace = workspace
                print(Feature_name)
                LAS[Feature_name] = file_name
                print(LAS[Feature_name])

            else:
                arcpy.AddMessage('Downloading data from %s'%url_save)
                workspace_zip = workspace + '\\' + title[:-4]
                urllib.request.urlretrieve(url_save, workspace + '\\' + title)
                # Unzip the data
                arcpy.AddMessage("Unzipping File...")
                if os.path.exists(workspace_zip):
                    print('File already exist...')
                else:
                    zip_ref = zipfile.ZipFile(workspace + '\\' + title)
                    title = title[:-4]
                    os.makedirs(workspace + '\\' + title)
                    workspace_zip = workspace + '\\' + title
                    print('workspace:', workspace)
                    print('workspace_zip', workspace_zip)
                    zip_ref.extractall(workspace_zip)
                    file_name = workspace_zip + '\\' + title + '.las'
                    if title[0] == 'N':
                        create_epsg(title[:-4], file_name)
                        arcpy.env.workspace = workspace
                    print(Feature_name)
                    LAS[Feature_name] = file_name
                    print(LAS[Feature_name])
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

            # Set moved to 0 indicating a new summit
            moved = 0

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
    arcpy.AddField_management(FCFile, 'source', 'TEXT')
    # cursor = arcpy.da.InsertCursor(FCFile, newFields)
    for entry in pointsToAddOverall:
        with arcpy.da.InsertCursor(FCFile, newFields) as cursor:
            row = [entry['NAME'], entry['FID'], entry['STATE'], entry['COUNTY'], entry['FEATURE'], entry['POINT'],
                    entry['DISTANCE'], entry['ELEV'], entry['SOURCE']]
            cursor.insertRow(row)
            arcpy.AddMessage(row)
    arcpy.AddMessage('Done!')

    return (FCFile)


def combine(infc):
    for feateures in arcpy.SearchCursor(infc):
        if feateures.getValue('FEATURE_NA') not in feature_name:
            feature_name.append(feateures.getValue('FEATURE_NA'))

    for summit in feature_name:
        print(summit)
        cursor = arcpy.SearchCursor(infc)
        for row in cursor:
            title = row.getValue('FEATURE_NA')
            pointFID = row.getValue('FEATURE_ID')
            pointState = row.getValue('STATE_ALPH')
            pointCounty = row.getValue('COUNTY_NAM')
            feature = row.getValue('FEATURE_CL')
            if row.getValue('source') == '1/3':
                elevation_10m = row.getValue('elevation')
                pointshape = row.getValue('Shape')
            print("FID:", pointFID)
            feature_next = cursor.next()
            while feature_next:
                print("While Loop", feature_next.getValue("FEATURE_NA"))
                if summit == feature_next.getValue('FEATURE_NA'):
                    print("In for loop,", pointFID, "==", feature_next.getValue('FEATURE_ID'))
                    if feature_next.getValue('source') == '1m':
                        elevation_1m = feature_next.getValue('elevation')
                        pointshape = feature_next.getValue('Shape')
                    if feature_next.getValue('source') == 'LPC':
                        elevation_LPC = feature_next.getValue('elevation')
                        pointshape = feature_next.getValue('Shape')
                feature_next = cursor.next()

        pointsToAddOverall.append({'NAME': title, 'FID': pointFID, 'STATE': pointState, 'COUNTY': pointCounty,
                                       'FEATURE': feature, 'POINT': pointshape, 'ELEV10M': elevation_10m, 'ELEV1M': elevation_1m,
                                       'ELEVLPC':  elevation_LPC})
        with arcpy.da.UpdateCursor(infc, "FEATURE_NA") as cursor:
            for row in cursor:
                if row[0] ==summit :
                    cursor.deleteRow()

    arcpy.AddMessage('Writing data to feature class....')
    # Add data to the new feature class
    FinalFields = ['FEATURE_NA', 'FEATURE_ID', 'STATE_ALPH', 'COUNTY_NAM','SHAPE@','FEATURE_CL', 'elev_10m', 'elev_1m', 'elev_LPC']
    FCFile = outpath + '\\' + outfile
    arcpy.DeleteField_management(FCFile,['STATW_NUME','COUNTY_NUM','PRIMARY_LA','PRIM_LONG_','PRIMARY_1','PRIM_LONG1','PRIM_LAT_D',
                                         'PRIM_LON_1','SOURCE_LAT','SOURCE_LON','SOURCE_L_1','SOURCE_L_2','SOURCE_L_3','SOURCE_L_4',
                                         'ELEV_IN_M', 'ELEV_IN_FT','MAP_NAME','DATE_CREAT','DATE_EDITE','Near_Pts', 'distance_m', 'source', 'sum_count', 'elevation'])
    arcpy.AddField_management(FCFile, 'elev_10m', 'FLOAT')
    arcpy.AddField_management(FCFile, 'elev_1m', 'FLOAT')
    arcpy.AddField_management(FCFile, 'elev_LPC', 'FLOAT')
    # cursor = arcpy.da.InsertCursor(FCFile, newFields)
    for entry in pointsToAddOverall:
        with arcpy.da.InsertCursor(FCFile, FinalFields) as cursor:
            row = [entry['NAME'], entry['FID'], entry['STATE'], entry['COUNTY'],entry['POINT'],entry['FEATURE'], entry['ELEV10M'], entry['ELEV1M'], entry['ELEVLPC']]
            cursor.insertRow(row)
            arcpy.AddMessage(row)
    arcpy.AddMessage('Done!')
'''-----------------------------Variables----------------------------------'''
# Variables for the user define functions

# Variable to keep track of how many pixels a summit moved
global moved
moved = 0
# temp[] is use to make the polygon to clip the shapefile
temp = []
temp2 = []
# bbox[] is use to download/get the srcLocation (TIF/IMG)
bbox = []
global merge_file
merge_file = []
# Bool used to exit the recursive loop
summitFound = False
# List to be filled with the summit found for each point analyzed
pointsToAddOverall = []
# final feature class
pointsfinals = []
# delete the feature that dont have raster data
global delete_feature
delete_feature = []
# Fields to retreive from "infc" and add to the result feature class
#fields = ['gaz_name', 'feature_id', 'state_alph', 'county_nam', 'gaz_featur', 'SHAPE@XY']
fields = ['FEATURE_NA', 'FEATURE_ID', 'STATE_ALPH', 'COUNTY_NAM', 'FEATURE_CL', 'SHAPE@']
#newFields = ['gaz_name', 'feature_id', 'state_alph', 'county_nam', 'gaz_featur', 'SHAPE@XY', 'distance']
newFields = ['FEATURE_NA', 'FEATURE_ID', 'STATE_ALPH', 'COUNTY_NAM', 'FEATURE_CL', 'SHAPE@', 'distance_m', 'elevation','Source']
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
# Original X and Y values
LAS ={}
global GNISx
global GNISy
GNISx = 0
GNISy = 0

# A checkbox to ask the user whether they want to download gnis_data
Checkbox_GNIS = arcpy.GetParameterAsText(0)
# Dir that would be saved for the GNIS shapefiles
Save_Dir_GNIS = arcpy.GetParameterAsText(1)
#Dir that will save the final result
Summit_Dir = arcpy.GetParameterAsText(5)
# Ask user for the infc, shapfile if they don't want to download the update
infc = arcpy.GetParameterAsText(3)
# Bounding box to creat the TIF files and polygon to clip the shapefile(infc) (xin,ymin,xmax,ymax) in NAEC
BBox_Value = arcpy.GetParameterAsText(2)
# Buffer for END testing, should be an INTcc
buffer = arcpy.GetParameter(4)
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
    arcpy.Clip_analysis(infc,polygon,infc+'_StudyArea.shp')
    infc = infc +'_StudyArea.shp'
else:
    arcpy.AddMessage("Selecting only summit points")
    out_name = infc + '_Summit.shp'
    arcpy.Select_analysis(infc, out_name,
                          '"FEATURE_CL" = \'Summit\'')
    infc = out_name

    arcpy.MinimumBoundingGeometry_management(infc, workspace+'\\poly.shp', geometry_type='RECTANGLE_BY_WIDTH')
    out_cor = arcpy.SpatialReference(6318)
    arcpy.AddGeometryAttributes_management(workspace+'\\poly.shp',Geometry_Properties='EXTENT',Coordinate_System=out_cor)
    with arcpy.da.SearchCursor(workspace+'\\poly.shp', ['EXT_MIN_X', 'EXT_MIN_Y','EXT_MAX_X','EXT_MAX_Y']) as cursor:
        for row in cursor:
            bbox.append({'NAME': 'GNIS_BBOX', 'POINT': (float(row[0]),float(row[3]))})
            bbox.append({'NAME': 'GNIS_BBOX', 'POINT': (float(row[2]), float(row[1]))})

for item in bbox:
    print(item)
# Create_feature will return the path for the clipped TIF/IMG
srcLocation = Create_feature('GNIS',bbox)
merge_file = []
xmin = float(arcpy.GetRasterProperties_management(srcLocation,'LEFT').getOutput(0))
ymax = float(arcpy.GetRasterProperties_management(srcLocation,'TOP').getOutput(0))
# Size of the cells in the raster image
CellSizeX = float(arcpy.GetRasterProperties_management(srcLocation, 'CELLSIZEX').getOutput(0))
CellSizeY = float(arcpy.GetRasterProperties_management(srcLocation, 'CELLSIZEY').getOutput(0))
# get the spatial_ref from the source location, should be North America Equidistant Coic
spatial_ref = arcpy.Describe(srcLocation).spatialReference

infc = workspace+'\\'+infc
arcpy.AddMessage('infc:%s'%infc)
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

print(srcLocation)
print(infc)
# Find the summit points for 1/3(10m) Raster, return the path to the shapefile(infc)
TenMeterSummit = FindSummit(srcLocation, infc, workspace,'TenMeterSummit')
print ('processing data...' )
# Merge the raster files for testing the one meter summit
spat_ref_NAEC = arcpy.SpatialReference('North America Equidistant Conic')
if len(merge_file) > 1:
    arcpy.AddMessage("Merge Raster...")
    arcpy.MosaicToNewRaster_management(merge_file,workspace,'\\MergedIMG_NAEC1m.tif',coordinate_system_for_the_raster=spat_ref_NAEC,
                                       pixel_type='32_BIT_FLOAT', number_of_bands=1,mosaic_method='LAST',mosaic_colormap_mode='FIRST')
    srcLocation = (workspace+'\\MergedIMG_NAEC1m.tif')
elif len(merge_file) == 0:
    arcpy.CopyFeatures_management(TenMeterSummit,Summit_Dir)
    arcpy.AddMessage("No one meter data...")
    arcpy.AddMessage("Done")
    quit()
else:
    print(merge_file)
    arcpy.MosaicToNewRaster_management(merge_file[0],workspace,'\\MergedIMG_NAEC1m.tif',coordinate_system_for_the_raster=spat_ref_NAEC,
                                       pixel_type='32_BIT_FLOAT', number_of_bands=1,mosaic_method='LAST',mosaic_colormap_mode='FIRST')
    srcLocation = (workspace + '\\MergedIMG_NAEC1m.tif')


merge_file = []
pointsToAddOverall = []
neighbors.clear()


# Size of the cells in the raster image
CellSizeX = float(arcpy.GetRasterProperties_management(srcLocation, 'CELLSIZEX').getOutput(0))
CellSizeY = float(arcpy.GetRasterProperties_management(srcLocation, 'CELLSIZEY').getOutput(0))
# Convert Raster to integer
OutRas = arcpy.sa.Con(Raster(srcLocation) > -999, 1, 0)
IntRas = arcpy.Int_3d(OutRas)
arcpy.AddMessage('Raster converted to integer raster')

# Convert the Raster to a polygon feature class
arcpy.RasterToPolygon_conversion(IntRas, tempSRC)
arcpy.AddMessage('Integer raster converted to a polygon feature class')

# Clip the feature class to the extent of the raster
arcpy.Clip_analysis(TenMeterSummit, tempSRC, clippedFC)
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

#code to delete the unwanted features
if len(delete_feature) >=1:
    if len(delete_feature) ==1:
        sql_exp = """{0} NOT IN '{1}'""".format(arcpy.AddFieldDelimiters(TenMeterSummit, 'FEATURE_NA'),delete_feature[0])
    else:
        sql_exp = """{0} NOT IN {1}""".format(arcpy.AddFieldDelimiters(TenMeterSummit, 'FEATURE_NA'),str(tuple(delete_feature)))
    arcpy.Select_analysis(TenMeterSummit, 'selected.shp', sql_exp)
    OneMeterSummit = FindSummit(srcLocation, 'selected.shp', workspace,'OneMeterSummit')
else:
    OneMeterSummit = FindSummit(srcLocation, TenMeterSummit, workspace, 'OneMeterSummit')


merge_file = []
pointsToAddOverall = []
neighbors.clear()


'''
Code to filter out the LPC points, select the point with in 2M radius and find the highest one
'''
arcpy.FeatureTo3DByAttribute_3d (OneMeterSummit,workspace+'\\OneMeterSummit3D.shp' , 'elevation')
in_feature = workspace+'\\OneMeterSummit3D.shp'

# Create a feature class for all the highest shape point
spatial_ref = arcpy.SpatialReference(102010)
arcpy.CreateFeatureclass_management(workspace, 'LPC', "POINT",in_feature, "DISABLED", "DISABLED", spatial_ref)

cursor = arcpy.SearchCursor(in_feature)


for row in cursor:
    title = row.getValue('FEATURE_NA')
    pointFID = row.getValue('FEATURE_ID')
    pointState = row.getValue('STATE_ALPH')
    pointCounty = row.getValue('COUNTY_NAM')
    feature = row.getValue('FEATURE_CL')
    Path_title = title.replace(" ", "_")
    sql = """{0} = '{1}'""".format(arcpy.AddFieldDelimiters(in_feature, 'FEATURE_NA'),title)
    arcpy.Select_analysis(in_feature, Path_title, sql)
    print(LAS[title])
    layer = arcpy.MakeLasDatasetLayer_management (LAS[title], out_layer='Ground', class_code=[2], return_values='Last Return')
    arcpy.LocateLasPointsByProximity_3d(layer, Path_title+'.shp', search_radius='2 Meter', count_field='Near_Pts',out_features= Path_title+'_Select', geometry='Point')
    in_file = Path_title+'_Select.shp'
    if arcpy.management.GetCount(in_file)[0] != "0":
        in_file_prj = Path_title+'SelectNAEC.shp'
        out_corsys = arcpy.SpatialReference(102010)
        arcpy.AddZInformation_3d(in_feature_class=in_file,out_property='Z')
        arcpy.Project_management(in_file, in_file_prj, out_corsys)
        arcpy.AddXY_management(in_file_prj)
        maxValue = 0
        summit_num = 1
        for entries in arcpy.SearchCursor(in_file_prj):
            print(entries.getValue('Z'))
            if maxValue == entries.getValue('Z'):
                summit_num += 1
            if maxValue < entries.getValue('Z'):
                maxValue = entries.getValue('Z')
                pointshape = entries.getValue('Shape')
                x = entries.getValue('POINT_X')
                print (x)
                y = entries.getValue('POINT_Y')
                print (y)
                print (GNISx, GNISy)
                a = abs(GNISx - x)
                b = abs(GNISy - y)
                distance = math.sqrt((pow(a, 2)) + (pow(b, 2)))


        pointsToAddOverall.append({'NAME': title, 'FID': pointFID, 'STATE': pointState, 'COUNTY': pointCounty,
                                                   'FEATURE': feature,'POINT':pointshape,'ELEV': maxValue,'SOURCE': 'LPC', 'Dist': distance,'Sum_count':summit_num})

arcpy.AddMessage('Writing data to feature class....')
# Add data to the new feature class
FinalFields = ['FEATURE_NA', 'FEATURE_ID', 'STATE_ALPH', 'COUNTY_NAM','SHAPE@','FEATURE_CL', 'elevation', 'source', 'distance_m', 'Sum_count']
FCFile = workspace + '\\' + 'LPC.shp'
arcpy.DeleteField_management(FCFile,['STATW_NUME','COUNTY_NUM','PRIMARY_LA','PRIM_LONG_','PRIMARY_1','PRIM_LONG1','PRIM_LAT_D',
                                     'PRIM_LON_1','SOURCE_LAT','SOURCE_LON','SOURCE_L_1','SOURCE_L_2','SOURCE_L_3','SOURCE_L_4',
                                     'ELEV_IN_M', 'ELEV_IN_FT','MAP_NAME','DATE_CREAT','DATE_EDITE','Near_Pts'])
arcpy.AddField_management(FCFile, 'source', 'TEXT')
arcpy.AddField_management(FCFile, 'sum_count', 'LONG')
# cursor = arcpy.da.InsertCursor(FCFile, newFields)
for entry in pointsToAddOverall:
    with arcpy.da.InsertCursor(FCFile, FinalFields) as cursor:
        row = [entry['NAME'], entry['FID'], entry['STATE'], entry['COUNTY'],entry['POINT'],entry['FEATURE'], entry['ELEV'], entry['SOURCE'], entry['Dist'], entry['Sum_count']]
        cursor.insertRow(row)
        arcpy.AddMessage(row)
arcpy.AddMessage('Done!')





input_file = [workspace+'\\TenMeterSummit.shp',workspace+'\\OneMeterSummit.shp',workspace+'\\LPC.shp']
arcpy.Merge_management (input_file, "combine.shp", field_mappings='JOIN')
arcpy.DeleteField_management(workspace + "\\combine.shp",['STATW_NUME','COUNTY_NUM','PRIMARY_LA','PRIM_LONG_','PRIMARY_1','PRIM_LONG1','PRIM_LAT_D',
                                     'PRIM_LON_1','SOURCE_LAT','SOURCE_LON','SOURCE_L_1','SOURCE_L_2','SOURCE_L_3','SOURCE_L_4',
                                     'ELEV_IN_M', 'ELEV_IN_FT','MAP_NAME','DATE_CREAT','DATE_EDITE'])
LAS.clear()

infc = workspace + "\\combine.shp"

arcpy.AddMessage(Summit_Dir)

spatial_ref = arcpy.SpatialReference(102010)
Summit_Dir = Summit_Dir.split('\\')
arcpy.AddMessage(Summit_Dir)
outpath = ''

for i in range(len(Summit_Dir)):
    if i != len(Summit_Dir)-1:
        outpath += Summit_Dir[i]
        if i != len(Summit_Dir)-2:
            outpath += "\\"
outfile = Summit_Dir[-1]

arcpy.CreateFeatureclass_management(outpath, outfile, "POINT",infc, "DISABLED", "DISABLED", spatial_ref)

feature_name = []
pointsToAddOverall = []

combine(infc)