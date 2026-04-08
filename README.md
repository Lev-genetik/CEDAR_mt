# Repository overview

Essential code and data of the Shagam, L. et al, "Multi-level aggregation analysis of microbiome composition and host gene expression reveals associations with systemic and local immunity", 2026 manuscript.

## Folder structure:

- `notebooks`: R Jupyter notebooks performing 9 gut microbiome-host gene expression comparisons at multiple feauture coarseness levels:
	- 02_mt.ipynb (4 comparisons)
	- 03_me_tc.ipynb (3 comparisons)
	- 04_en_t.ipynb (2 comparisons)

- `Data`: Input data for the notebooks. 
- `code`: Auxilliary functions and input data for the notebooks.
- `liege_mt.yaml`: Recipe for creating the Conda environment that cen be used to execute the notebooks (via VS Code, etc).

---

## Instructions:

All R notebooks are available with prerendered output. The main results can be reproduced by installing a Conda environment using the yaml file and notebooks and running the notebooks.

To set up the environment for running the notebooks, you will need:
- Conda
- Visual Studio Code (or similar IDE)

Create the Conda environment using the .yaml file that contains the required packages with their versions:
```bash
conda env create -f liege_mt.yaml
conda activate liege_mt
```

Populate the Data/CEDAR1_transcriptome folder with the input data files by downloading the archive here and unzipping it to that folder: https://drive.google.com/file/d/1KmAFmGMmkhi4AI4QGU83qV3jb0FY00UL/view?usp=sharing

The notebooks can be executed via VS Code, etc. using an R Jupyter kernel.

