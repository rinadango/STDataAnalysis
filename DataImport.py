# -*- coding: utf-8 -*-
"""
@author: Paulina Frolovaite

Two main methods: spatial_data_importing_complexdata() and spatial_data_importing_simpledata
"""

import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import anndata as ad
import numpy as np

import json
from PIL import Image
from numpy import asarray

# Paths to coordinate, scale factor and image files
def _path_spatial_folder(spatial_folder_path):
    file_list = ['tissue_positions_list.csv', 'scalefactors_json.json', 'tissue_hires_image.png', 'tissue_lowres_image.png']
    new_pathS = []
    for item in file_list:
        paths = f"{spatial_folder_path}/{item}"
        new_path = paths.replace('\\','/')
        new_pathS.append(new_path)
    return new_pathS

# IMAGE 
# GETTING IMAGE DATA IN NUMPY ARRAY
def _image_import(new_pathS):

    # hires
    # load the image
    II_image = Image.open(new_pathS[2])
    # convert image to numpy array
    image_data = asarray(II_image)

    #lowres
    low_II_image = Image.open(new_pathS[3])
    low_image_data = asarray(low_II_image)

    return image_data, low_image_data

# Getting the spatial factors
def _spatial_json(new_pathS):
    
    ### Opening JSON file
    f = open(new_pathS[1])
  
    # returns JSON object as 
    # a dictionary
    data_scale = json.load(f)

    # Closing file
    f.close()

    return data_scale

# Workaround the ANNDATA object
def _anndata_object_workaround(library_id, data_spatial, coord_df, data_scale, image_data, spatial_folder_path, image_type = None, low_image_data = None, both = False):

    # Key for ANNDATA .uns and .obsm
    spatial_key = "spatial"

    # Create dataset to add to ANNDATA.obs
    csv_3colms = coord_df.iloc[:, 0:4]
    csv_3colms = csv_3colms.set_index(csv_3colms.columns[0])
    csv_3colms.index.name = None
    csv_3colms.columns = ['in_tissue','array_row','array_col']
    
    # Get coordinates
    csv_COORDS = coord_df.iloc[:, [5, 4]].values

    # Add filtered out coordinates to ANNDATA.obsm
    data_spatial.obsm[spatial_key] = csv_COORDS
    
    # Creating ANNDATA.uns
    data_spatial.uns[spatial_key] = {library_id: {}}
    data_spatial.uns[spatial_key][library_id]["images"] = {}

    # Add scaling factors to ANNDATA.uns
    data_spatial.uns[spatial_key][library_id]["scalefactors"] = data_scale
    
    # ADDING SOURCE IMAGE PATH TO ANNDATA.UNS
    data_spatial.uns[spatial_key][library_id]['metadata'] = {'source_image_path': spatial_folder_path}
    
    # Want hires or lowres
    if image_type == "lowres":
        # ADDING IMAGE DATA (NUMPY ARRAY) TO ANNDATA.UNS
        data_spatial.uns[spatial_key][library_id]["images"] = {"lowres": image_data}
    elif image_type == "hires":
        data_spatial.uns[spatial_key][library_id]["images"] = {"hires": image_data}
    else:
        pass

    if both == True:
        data_spatial.uns[spatial_key][library_id]["images"] = {"hires": image_data, "lowres": low_image_data}
    else:
        pass
    
    # ADD ['in_tissue','array_row','array_col'] to ANNDATA.OBS
    data_spatial.obs = csv_3colms

    return data_spatial

"""
USE `spatial_data_importing_complexdata()` only for data that has `mtx.gz` and `tsv.gz` files.

Args:
    library_id: ID of the library; INPUT `string`
    data_path: file path to the `matrix.mtx.gz`, `barcodes.tsv.gz` and `features.tsv.gz` files folder. INPUT is a `file path`
    spatial_folder_path: path to the `spatial` folder; INPUT `file path`
    imag_type: specify either as `lowres` or `hires`
    both: True or False
Returns:
    data_spatial: ANNDATA object
    
Paths can have `\\` or `/` depending on the OS used
For WIN: it is important to change the `\` to `/` or `\\`, or add `r` at the begining of the `string`. RECOMMENDATION: to add the `r` at the begining of the `string`

Takes either both high- and low-resolution images, or either
"""

def spatial_data_importing_complexdata(library_id, data_path, spatial_folder_path, image_type = None, both = False):
    
    # Import matrix, barcodes, features
    data_spatial = sc.read_10x_mtx(data_path)
    
    data_spatial.var_names_make_unique()
    
    # Get spatial folder path
    new_pathS = _path_spatial_folder(spatial_folder_path)
    
    ### Import coordinates file
    coord_df = pd.read_csv(new_pathS[0], sep=',',header=None)
    
    # Read JSON
    data_scale = _spatial_json(new_pathS)
    
    # HIRES and LOWRES Import image data
    hi_image_data, low_image_data = _image_import(new_pathS)
    
    # Add everything into ANNDATA object
    if image_type == 'hires' and both == False:
        data_spatial = _anndata_object_workaround(library_id, data_spatial, coord_df,
        data_scale, hi_image_data, spatial_folder_path, image_type = 'hires', both = False)
    elif image_type == 'lowres' and both == False:
        data_spatial = _anndata_object_workaround(library_id, data_spatial, coord_df,
        data_scale, low_image_data, spatial_folder_path, image_type = 'lowres', both = False)
    else:
        pass

    if both == True and image_type == None:
        data_spatial = _anndata_object_workaround(library_id, data_spatial, coord_df,
        data_scale, hi_image_data, spatial_folder_path, low_image_data = low_image_data, both = True)
    elif both == False and image_type == None:
        raise ValueError("Parameter 'both' must not be false without specifying image type. Either specify 'both' as True or specify image type.")

    return data_spatial


def spatial_data_importing_simpledata(library_id, matrix, barcodes_filtered, feature_file, spatial_folder_path, image_type = None, both = False):
    
    """
    USE `spatial_data_importing_simpledata()` only for data that IS a simple tab, space or comma delimited file. Basically for files that are `tsv`, `csv`, `txt` and so on.

    Args:
        library_id: ID of the library; INPUT `string`
        matrix: count data file taken from the `filtered_count_matrices` folder. INPUT is`file path`
        feature_file: features file; INPUT: `file path`; NOT zipped
        barcodes_filtered: `barcodes` file taken from the `filtered_count_matrices` folder. INPUT is a `file path`. NOT zipped
        spatial_folder_path: path to the `spatial` folder. INPUT is a `file path`
        imag_type: specify either as `lowres` or `hires`
        both: True or False
    Returns:
        data_spatial: ANNDATA object

    Paths can have `\\` or `/` depending on the OS used
    
    For WIN: it is important to change the `\` to `/` or `\\`, or add `r` at the begining of the `string`. RECOMMENDATION: to add the `r` at the begining of the `string`
    """

    ### Import matrix
    data_spatial = sc.read(matrix, cache=True)
    # Fix shape
    data_spatial = data_spatial.transpose()
    
    # Get spatial folder path
    new_pathS = _path_spatial_folder(spatial_folder_path)
    
    ### Importing
    # Import coordinates file
    coord_df = pd.read_csv(new_pathS[0], sep=',',header=None)
    # Import barcodes file
    barcodes_filtered_df = pd.read_csv(barcodes_filtered, header = None)
    # Import features file
    with open(feature_file, 'r') as f:
        featureS = [line.strip() for line in f]
        f.close()
        
    # Creating ANNDATA.var
    data_spatial.var_names = featureS

    # Filter coordinates file based on barcodes file
    filtered_out_csv = pd.merge(coord_df, barcodes_filtered_df, how='inner')

    # Read JSON
    data_scale = _spatial_json(new_pathS)
  
    # Import image data
    hi_image_data, low_image_data = _image_import(new_pathS)

   # Add everything into ANNDATA object
    if image_type == 'hires' and both == False:
        data_spatial = _anndata_object_workaround(library_id, data_spatial, filtered_out_csv,
        data_scale, hi_image_data, spatial_folder_path, image_type = 'hires', both = False)
    elif image_type == 'lowres' and both == False:
        data_spatial = _anndata_object_workaround(library_id, data_spatial, filtered_out_csv,
        data_scale,  low_image_data, spatial_folder_path, image_type = 'lowres', both = False)

    if both == True and image_type == None:
        data_spatial = _anndata_object_workaround(library_id, data_spatial, filtered_out_csv,
        data_scale, hi_image_data, spatial_folder_path, image_type = None, both = True , low_image_data = low_image_data)
    elif both == False and image_type == None:
        raise ValueError("Parameter 'both' must not be false without specifying image type. Either specify 'both' as True or specify image type.")

    return data_spatial

