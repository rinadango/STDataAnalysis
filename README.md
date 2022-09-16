# Usage

<b>Important</b>: must have scanpy and anndata installed for `DataImport.py`<p>
Add the folder where you will use it as import for other code<br>
Import main module as: `from STDataAnalysis import <.py file name>`<p>
Read comments of each `.py` for more instructions

## DataImport.py

USE `spatial_data_importing_raw` only for data files that are `mtx.gz`, `tsv.gz` or other more complex files. <br>
USE `spatial_data_importing_filtered` only for data files are simple tab, space or comma delimited files. Basically for files that are `tsv`, `csv`, `txt` and so on.

## FilterBYName.py

USE `filter_by_type_gene` for filtering genes by their type