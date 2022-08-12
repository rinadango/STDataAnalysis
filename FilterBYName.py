# -*- coding: utf-8 -*-
"""
@author: Paulina Frolovaite

`gene_starts` - identifier with what the gene name starts, for eg.: MT- for mitochondrials genes; INPUT: a string or a list of strings
`anndata_object` - the imported data set in an anndata object having features
"""

def filter_by_type_gene(gene_starts, anndata_object):

    anndata_object = [names for names in anndata_object.var_names]

    type_of_gene = type(gene_starts) # check type of variable

    if type_of_gene == str: # if just one indentifier
        for gene_id in anndata_object:
            if gene_id.startswith(gene_starts):
                print(gene_id)

    elif type_of_gene == list: # if several
        for item in gene_starts:
            for gene_id in anndata_object:
                if gene_id.startswith(item):
                    print(gene_id)