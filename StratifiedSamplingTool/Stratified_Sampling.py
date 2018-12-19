import arcpy

'''
This script will generate a stratified sample for testing
Arthur Chan
'''

In_shapefile = arcpy.GetParameterAsText(0)
#if checkbox == false, create random point, locations from user provid point file otherwise
CheckBox_pointfile = arcpy.GetParameterAsText(1)
Point_file = arcpy.GetParameterAsText(2)
Accuracy = arcpy.GetParameter(3)
Error_Margin = arcpy.GetParameter(4)
Save_file = arcpy.GetParameterAsText(5)
Total_point = arcpy.GetParameter(6)
Workspace = arcpy.GetParameterAsText(7)

arcpy.env.workspace = Workspace
out_cor = arcpy.SpatialReference(102008)
global point
point = Total_point

def calculate_point_per_poly(shape_file):
    total_area = 0
    arcpy.AddMessage('Generating area for features...')
    arcpy.AddField_management(shape_file,field_name='Temp_Area',field_type='DOUBLE')
    arcpy.CalculateGeometryAttributes_management(shape_file,geometry_property=[['Temp_Area','AREA_GEODESIC']]
                                                 ,area_unit='SQUARE_METERS', coordinate_system=out_cor)

    #get the total area
    arcpy.AddMessage('Computing the total area...')
    for row in arcpy.da.SearchCursor(shape_file,'Temp_Area'):
        total_area += row[0]
        #print('Temp_Area',row[0])

    #calculate the sampling size via total are and distribute to each feature according to area_percentage
    arcpy.AddField_management(shape_file,field_name='SampleSize',field_type='SHORT')
    arcpy.AddMessage('Computing the sampling size for each feature...')
    with arcpy.da.UpdateCursor(shape_file,['Temp_Area','SampleSize']) as cursor:

        for row in cursor:
            area_percent = row[0]/total_area
            if CheckBox_pointfile == 'false':
                row[1]=int(round(int(point)*area_percent))
            elif CheckBox_pointfile == 'true':
                n = (1.96 * 1.96 * int(Accuracy) * (100 - int(Accuracy))) / (int(Error_Margin) * int(Error_Margin))
                Total_point = arcpy.GetCount_management(Point_file)
                total_point = n / (1 + ((n - 1) / int(Total_point[0])))
                row[1] = int(round(total_point * area_percent))
            #print('Sample Size:',row[1])
            cursor.updateRow(row)

    #print('total area:',total_area)

#Generate sampling points for each feature
def generate_point(shape_file):
    merge = []
    for row in arcpy.da.SearchCursor(shape_file, ['FID','SampleSize']):
        sql = """{0} = {1}""".format('FID', row[0])
        out_name = 'FID'+str(row[0])
        constrain_feat = arcpy.Select_analysis(in_features=shape_file,out_feature_class=out_name,where_clause=sql)
        if CheckBox_pointfile == 'false':
            arcpy.AddMessage('Creating random points for features...')
            out_name += '_testpoints'
            arcpy.CreateRandomPoints_management(out_path=Workspace, out_name=out_name
                                                ,constraining_feature_class=constrain_feat
                                                ,number_of_points_or_field=row[1],minimum_allowed_distance='1 Meter'
                                                ,create_multipoint_output='POINT')
            merge.append(out_name+'.shp')
        else:
            #select points within the polygon
            arcpy.AddMessage('Selecting random points for features...')
            in_feat= arcpy.SelectLayerByLocation_management(in_layer=Point_file,overlap_type='INTERSECT'
                                                                ,select_features=out_name+'.shp', selection_type='NEW_SELECTION')
            out_name += '_testpoints'
            #print(in_feat)
            arcpy.SubsetFeatures_ga(in_features=in_feat,out_training_feature_class=out_name
                                    ,size_of_training_dataset=row[1],subset_size_units='ABSOLUTE_VALUE')
            merge.append(out_name+'.shp')
    #print(merge)
    arcpy.AddMessage('Merging files...')
    arcpy.Merge_management(merge, Save_file)

    arcpy.AddMessage('Points generated')

''' ---------- Main ---------- '''
arcpy.env.overwriteOutput = True

calculate_point_per_poly(In_shapefile)
generate_point(In_shapefile)
arcpy.DeleteField_management(In_shapefile,['SampleSize','Temp_Area'])
arcpy.AddMessage('Done')