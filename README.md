# python_ps

[![DOI](https://zenodo.org/badge/76561413.svg)](https://zenodo.org/badge/latestdoi/76561413)

Welcome to the repository of PlanetServer's Python API. This respository contains files associated with the project. A small description of each file is given below:

- metadata.txt : text file containing the different band names and their corresponding wavelengths

- bands_table.ipynb : Jupyter Notebook containing the python code to extract the band names and wavelengths from metadata.txt 

- bands_table.fits : FITS file (astropy standards) containing a table with the extracted band names and wavelengths

- summary_products.ipynb : Jupyter Notebook integrating CRISM summary products and URL generator (WCPS queries) for different RGB band math combinations

- script_summary_products.ipynb : summary_products.ipynb merged in a single cell --> equivalent to a python script
