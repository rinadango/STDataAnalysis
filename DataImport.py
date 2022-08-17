# -*- coding: utf-8 -*-
"""
@author: Paulina Frolovaite

Two main methods: spatial_data_importing_raw() and spatial_data_importing_filtered
"""

import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import anndata as ad
import squidpy as sq
import numpy as np

import json
from PIL import Image
from numpy import asarray

# Paths to coordinate, scale factor and image files
def _path_spatial_folder(spatial_folder_path):
    file_list = ['tissue_positions_list.csv', 'scalefactors_json.json', 'tissue_hires_image.png']
    new_pathS = []
    for item in file_list:
        paths = f"{spatial_folder_path}/{item}"
        new_path = paths.replace('\\','/')
        new_pathS.append(new_path)
    return new_pathS

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

# IMAGE 
# GETTING IMAGE DATA IN NUMPY ARRAY
def _image_import(new_pathS):

    # load the image
    II_image = Image.open(new_pathS[2])
    # convert image to numpy array
    image_data = asarray(II_image)

    return image_data

# Workaround the ANNDATA object
def _anndata_object_workaround(library_id, data_spatial, coord_df, data_scale, image_data, spatial_folder_path):

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
    
    # ADDING IMAGE DATA (NUMPY ARRAY) TO ANNDATA.UNS
    data_spatial.uns[spatial_key][library_id]["images"] = {"hires": image_data}
    
    # ADD ['in_tissue','array_row','array_col'] to ANNDATA.OBS
    data_spatial.obs = csv_3colms

    return data_spatial

"""
`library_id` - ID of the library; INPUT `string`
`data_path` - file path to the `matrix.mtx.gz`, `barcodes.tsv.gz` and `features.tsv.gz` files folder; INPUT `file path`
`spatial_folder_path` - path to the `spatial` folder; INPUT `file path`
    
#### Paths can have `\\` or `/` depending on the OS used
For WIN: it is important to change the `\` to `/` or `\\`, or add `r` at the begining of the `string`. RECOMMENDATION: to add the `r` at the begining of the `string`
    
### USE `spatial_data_importing_raw()` only for RAW data

### Only takes high resolution image for now
"""

def spatial_data_importing_raw(library_id, data_path, spatial_folder_path):
    
    # Import matrix, barcodes, features
    data_spatial = sc.read_10x_mtx(data_path)
    
    data_spatial.var_names_make_unique()
    
    # Get spatial folder path
    new_pathS = _path_spatial_folder(spatial_folder_path)
    
    ### Import coordinates file
    coord_df = pd.read_csv(new_pathS[0], sep=',',header=None)
    
    # Read JSON
    data_scale = _spatial_json(new_pathS)
    
    # Import image data
    image_data = _image_import(new_pathS)
    
    # Add everything into ANNDATA object
    data_spatial = _anndata_object_workaround(library_id, data_spatial, coord_df,
    data_scale, image_data, spatial_folder_path)

    return data_spatial


"""
`library_id` - ID of the library; INPUT `string`
`matrix` - count data file taken from the `filtered_count_matrices` folder; INPUT: `file path`; available extensions `{'tab', 'soft.gz', 'loom', 'mtx', 'anndata', 'h5', 'txt', 'data', 'h5ad', 'csv', 'tsv', 'mtx.gz', 'xlsx'}`
`feature_file` - features file; INPUT: `file path`; NOT zipped
`barcodes_filtered`- `barcodes` file taken from the `filtered_count_matrices` folder; INPUT: `file path`; NOT zipped
`spatial_folder_path` - path to the `spatial` folder; INPUT `file path`

#### Paths can have `\\` or `/` depending on the OS used
For WIN: it is important to change the `\` to `/` or `\\`, or add `r` at the begining of the `string`. RECOMMENDATION: to add the `r` at the begining of the `string`
    
### USE `spatial_data_importing_filtered()` only for FILTERED data

### Only takes in high resolution image for now
"""

def spatial_data_importing_filtered(library_id, matrix, barcodes_filtered, feature_file, spatial_folder_path):
 
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
    image_data = _image_import(new_pathS)
   
    # Add everything into ANNDATA object
    data_spatial = _anndata_object_workaround(library_id, data_spatial, filtered_out_csv,
    data_scale, image_data, spatial_folder_path)
    
    return data_spatial