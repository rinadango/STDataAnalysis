# -*- coding: utf-8 -*-
"""
@author: Paulina Frolovaite

`gene_starts` - identifier with what the gene name starts, for eg.: MT- for mitochondrials genes; INPUT: a string or a list of strings
`anndata_object` - the imported data set in an anndata object with features
"""

def filter_by_type_gene(gene_starts, anndata_object):

    anndata_object = [names for names in anndata_object.var_names]

    if type(gene_starts) == str: # if just one indentifier
        for gene_id in anndata_object:
            if gene_id.startswith(gene_starts):
                print(gene_id)

    elif type(gene_starts) == list: # if several
        for item in gene_starts:
            for gene_id in anndata_object:
                if gene_id.startswith(item):
                    print(gene_id)