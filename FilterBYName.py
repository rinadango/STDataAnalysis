# -*- coding: utf-8 -*-
"""
@author: Paulina Frolovaite

`gene_starts` - identifier with what the gene name starts, for eg.: MT- for mitochondrials genes; INPUT: a string
`anndata_object` - the imported data set in an anndata object with features
`calc_percentage` - optional argument for calculating ratio of certain gene ids
`amount`- optional argument for calculating amount of certain genes
"""

def _var_names(anndata_object):
    anndata_object = [names for names in anndata_object.var_names]
    return anndata_object

def _percentage_calc(all_ids, filtered_ids):

    filtered_ids = len(filtered_ids)
    all_ids = len(all_ids)

    return (filtered_ids/all_ids) * 100

def filter_by_type_gene(gene_starts: str, anndata_object, calc_percentage = False, amount = False):

    anndata_object = _var_names(anndata_object)

    filtered_ids = []
    
    for gene_id in anndata_object:
        if gene_id.startswith(gene_starts):
            filtered_ids.append(gene_id)

    if calc_percentage == True:
        print(f"Ratio of {gene_starts} genes in dataset (%):", _percentage_calc(anndata_object, filtered_ids))

    if amount == True:
        print(f"Number {gene_starts} genes in dataset:", len(filtered_ids))

    return filtered_ids