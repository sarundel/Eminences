import os.path
import urllib
import os
import requests
import arcpy
import urllib.request
import zipfile
import glob
import heapq


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

#Create feature and pass a bbox to Download data from tnm api
def Create_feature(TestData, list = []):
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
            bbox_create.append(row[0])
            bbox_create.append(row[1])
    # swap the position of ymax and ymin to fit the API
    temp_swap = bbox_create[1]
    bbox_create[1] = bbox_create[3]
    bbox_create[3] = temp_swap

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
            arcpy.AddMessage("Merge Raster...")
            arcpy.MosaicToNewRaster_management(merge_file, workspace, '\\MergedIMG_NAEC.tif',
                                               coordinate_system_for_the_raster=spat_ref_NAEC, pixel_type='32_BIT_FLOAT',
                                               number_of_bands=1, mosaic_method='LAST', mosaic_colormap_mode='FIRST')
            return (workspace + '\\MergedIMG_NAEC.tif')
        elif len(merge_file) == 0:
            print ('No Raster IMG found...')
        else:
            arcpy.MosaicToNewRaster_management(merge_file[0], workspace, '\\MergedIMG_NAEC.tif',
                                               coordinate_system_for_the_raster=spat_ref_NAEC, pixel_type='32_BIT_FLOAT',
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
    else:
        for i in range(0,down_len):
            url_save = data['items'][i]['downloadURL']
            print("Downloading info from requested web...")
            title = url_save.split('/')
            title = title[-1]
            title = title.replace(" ", "_")

            if os.path.exists(workspace + '\\' + title):
                workspace_zip = workspace + '\\' + title[:-4]
                print(workspace_zip)
                print (title + ': Already exists')
                if (workspace_zip + '\\*.img') in merge_file:
                    print(title + 'Already in list')
                else:
                    if len(merge_file) == 0:
                        merge_file.append(glob.glob(workspace_zip + '\\*.img'))
                    else:
                        merge_file.insert(0, glob.glob(workspace_zip + '\\*.img'))
            else:
                print ('Downloading data from:', url_save)
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
                    zip_ref.extractall(workspace_zip)
                    if len(merge_file) == 0:
                        merge_file.append(glob.glob(workspace_zip+'\\*.img'))
                    else:
                        merge_file.insert(0, glob.glob(workspace_zip+'\\*.img'))
                    zip_ref.close()



'----------------------------------Variables----------------------------------'
infc = arcpy.GetParameterAsText(1)          #GNIS.shp file
save_dir = arcpy.GetParameterAsText(2)      #dir for the new feature class
save_name = arcpy.GetParameterAsText(3)     #name for the new feature class
workspace = arcpy.GetParameterAsText(4)
arcpy.env.workspace = workspace

bbox = []                                   #To store the x,y-coordinate for download DEM
merge_file = []
pointsToAdd = []
fields = ['FEATURE_NA', 'FEATURE_ID', 'STATE_ALPH', 'COUNTY_NAM', 'FEATURE_CL', 'SHAPE@']
newFields = ['FEATURE_NA', 'FEATURE_ID', 'STATE_ALPH', 'COUNTY_NAM', 'FEATURE_CL', 'SHAPE@', 'elevation']
# Priority Queue used to iterate through the neighbors based on their values
neighbors = PQueue()

'----------------------------------Main----------------------------------'
arcpy.env.overwriteOutput = True

arcpy.AddMessage("Selecting only summit points")
arcpy.Select_analysis(infc, infc[:-4]+'_Summit.shp','"FEATURE_CL" = \'Summit\'')
infc = infc[:-4]+'_Summit.shp'

# get a bound box value from the shp file by creating a polygon
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
# "GNIS" will return 10m raster
srcLocation = Create_feature('GNIS',bbox)
xmin = float(arcpy.GetRasterProperties_management(srcLocation,'LEFT').getOutput(0))
ymax = float(arcpy.GetRasterProperties_management(srcLocation,'TOP').getOutput(0))

# get the spatial reference from the DEM
spatial_ref = arcpy.Describe(srcLocation).spatialReference

# Get value fir GNIS point base on the 10m raster
for row in arcpy.da.SearchCursor(infc, fields, spatial_reference=spatial_ref):
    pointName = row[0]
    pointFID = row[1]
    pointState = row[2]
    pointCounty = row[3]
    pointShape = row[5]
    x, y = pointShape.firstPoint.X, pointShape.firstPoint.Y
    feature = row[4]
    point = str(x) + ' ' + str(y)
    pointsToAdd.append({'NAME': pointName, 'FID': pointFID, 'STATE': pointState, 'COUNTY': pointCounty,
                        'FEATURE': feature, 'POINT': pointShape,
                        'ELEV': float(arcpy.GetCellValue_management(srcLocation, point).getOutput(0))})

arcpy.CreateFeatureclass_management(save_dir, save_name, "POINT", infc, "DISABLED", "DISABLED", spatial_ref)
arcpy.AddMessage('New feature class created')

arcpy.AddMessage('Writing data to feature class....')
# Add data to the new feature class
FCFile = save_dir + '\\' + save_name
if not FCFile.endswith('.shp'):
    FCFile = FCFile + '.shp'
arcpy.AddField_management(FCFile, 'Elevation', 'FLOAT')
# cursor = arcpy.da.InsertCursor(FCFile, newFields)
for entry in pointsToAdd:
    with arcpy.da.InsertCursor(FCFile, newFields) as cursor:
        row = [entry['NAME'], entry['FID'], entry['STATE'], entry['COUNTY'], entry['FEATURE'], entry['POINT'], entry['ELEV']]
        cursor.insertRow(row)
    arcpy.AddMessage(row)
arcpy.AddMessage('Done!')