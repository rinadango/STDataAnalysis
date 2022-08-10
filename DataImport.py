# -*- coding: utf-8 -*-
"""
@author: Paulina Frolovaite
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

"""
`library_id` - ID of the library; INPUT `string`
`data_path` - file path to the `matrix.mtx.gz`, `barcodes.tsv.gz` and `features.tsv.gz` files folder; INPUT `file path`
`spatial_folder_path` - path to the `spatial` folder; INPUT `file path`
    
#### Paths can have `\\` or `/` depending on the OS used
For WIN: it is important to change the `\` to `/` or `\\`, or add `r` at the begining of the `string`. RECOMMENDATION: to add the `r` at the begining of the `string`
    
### USE `spatial_data_importing_raw()` only for RAW data
"""

def spatial_data_importing_raw(library_id, data_path, spatial_folder_path):
    
    # Import matrix, barcodes, features
    data_spatial = sc.read_10x_mtx(data_path)
    
    data_spatial.var_names_make_unique()
    
    # Paths to coordinate, scale factor and image files
    file_list = ['tissue_positions_list.csv', 'scalefactors_json.json', 'tissue_hires_image.png']
    new_pathS = []
    for item in file_list:
        paths = f"{spatial_folder_path}\\{item}"
        new_path = paths.replace('\\','/')
        new_pathS.append(new_path)
    
    ### Import coordinates file
    coord_df = pd.read_csv(new_pathS[0], sep=',',header=None)
    
    # Create data set to add to ANNDATA.obs
    csv_3colms = coord_df.iloc[:, 0:4]
    csv_3colms = csv_3colms.set_index(csv_3colms.columns[0])
    csv_3colms.index.name = None
    csv_3colms.columns = ['in_tissue','array_row','array_col']
    
    # Get coordinates
    csv_COORDS = coord_df.iloc[:, [5, 4]].values
    
    # Key for ANNDATA .uns and .obsm
    spatial_key = "spatial"

    # Add filtered out coordinates to ANNDATA.obsm
    data_spatial.obsm[spatial_key] = csv_COORDS

    # Creating ANNDATA.uns
    data_spatial.uns[spatial_key] = {library_id: {}}
    data_spatial.uns[spatial_key][library_id]["images"] = {}
    
    ### Opening JSON file
    f = open(new_pathS[1])
  
    # returns JSON object as 
    # a dictionary
    data_scale = json.load(f)
    data_scale

    # Closing file
    f.close()
    
    # Add scaling factors to ANNDATA.uns
    data_spatial.uns[spatial_key][library_id]["scalefactors"] = data_scale
    
    # ADDING SOURCE IMAGE PATH TO ANNDATA.UNS
    data_spatial.uns[spatial_key][library_id]['metadata'] = {'source_image_path': spatial_folder_path}
    
    ### IMAGE 
    ## GETTING IMAGE DATA IN NUMPY ARRAY
    # load the image
    II_image = Image.open(new_pathS[2])
    # convert image to numpy array
    image_data = asarray(II_image)
    
    # ADDING IMAGE DATA (NUMPY ARRAY) TO ANNDATA.UNS
    data_spatial.uns[spatial_key][library_id]["images"] = {"hires": image_data}
    
    # ADD ['in_tissue','array_row','array_col'] to ANNDATA.OBS
    data_spatial.obs = csv_3colms
    
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
"""

def spatial_data_importing_filtered(library_id, matrix, barcodes_filtered, feature_file, spatial_folder_path):
 
    ### Import matrix
    data_spatial = sc.read(matrix, cache=True)
    # Fix shape
    data_spatial = data_spatial.transpose()
    
    # Paths to coordinate, scale factor and image files
    file_list = ['tissue_positions_list.csv', 'scalefactors_json.json', 'tissue_hires_image.png']
    new_pathS = []
    for item in file_list:
        paths = f"{spatial_folder_path}\\{item}"
        new_path = paths.replace('\\','/')
        new_pathS.append(new_path)
    
    ### Importing
    # Import coordinates file
    coord_df = pd.read_csv(new_pathS[0], sep=',',header=None)
    # Import barcodes file
    barcodes_filtered_df = pd.read_csv(barcodes_filtered, header = None)
    # Import features file
    with open(feature_file, 'r') as f:
        featureS = [line.strip() for line in f]
        f.close()
        
    # Filter coordinates file based on barcodes file
    filtered_out_csv = pd.merge(coord_df, barcodes_filtered_df, how='inner')
    
    # Create data set to add to ANNDATA.obs
    filtered_out_csv_3colms = filtered_out_csv.iloc[:, 0:4]
    filtered_out_csv_3colms = filtered_out_csv_3colms.set_index(filtered_out_csv_3colms.columns[0])
    filtered_out_csv_3colms.index.name = None
    filtered_out_csv_3colms.columns = ['in_tissue','array_row','array_col']
    
    # Get filtered out coordinates
    filtered_out_csv_COORDS = filtered_out_csv.iloc[:, [5, 4]].values
    
    # Creating ANNDATA.var
    data_spatial.var_names = featureS
    
    # Key for ANNDATA .uns and .obsm
    spatial_key = "spatial"
    
    # Add filtered out coordinates to ANNDATA.obsm
    data_spatial.obsm[spatial_key] = filtered_out_csv_COORDS
    
    # Creating ANNDATA.uns
    data_spatial.uns[spatial_key] = {library_id: {}}
    data_spatial.uns[spatial_key][library_id]["images"] = {}
    
    ### Opening JSON file
    f = open(new_pathS[1])
  
    # returns JSON object as 
    # a dictionary
    data_scale = json.load(f)
    data_scale

    # Closing file
    f.close()
    
    # Add scaling factors to ANNDATA.uns
    data_spatial.uns[spatial_key][library_id]["scalefactors"] = data_scale
    
    # ADDING SOURCE IMAGE PATH TO ANNDATA.UNS
    data_spatial.uns[spatial_key][library_id]['metadata'] = {'source_image_path': spatial_folder_path}
    
    ### IMAGE 
    ## GETTING IMAGE DATA IN NUMPY ARRAY
    # load the image
    II_image = Image.open(new_pathS[2])
    # convert image to numpy array
    image_data = asarray(II_image)
    
    # ADDING IMAGE DATA (NUMPY ARRAY) TO ANNDATA.UNS
    data_spatial.uns[spatial_key][library_id]["images"] = {"hires": image_data}
    
    # ADD ['in_tissue','array_row','array_col'] to ANNDATA.OBS
    data_spatial.obs = filtered_out_csv_3colms
    
    return data_spatial