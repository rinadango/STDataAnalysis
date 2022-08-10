# -*- coding: utf-8 -*-
"""
@author: Paulina Frolovaite

`gene_starts` - identifier with what the gene name starts, for eg.: MT- for mitochondrials genes; INPUT: a string or a list of strings
`variable_names` - a list of gene names/ids; INPUT: list
"""

def filter_by_type_gene(gene_starts, variable_names: list):

    type_of_gene = type(gene_starts) # check type of variable

    if type_of_gene == str: # if just one indentifier
        for gene_id in variable_names:
            if gene_id.startswith(gene_starts):
                print(gene_id)

    elif type_of_gene == list: # if several
        for item in gene_starts:
            for gene_id in variable_names:
                if gene_id.startswith(item):
                    print(gene_id)