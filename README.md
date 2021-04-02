# CAPD_CHIL2021 dataset preparation scripts for ["RNA Alternative Splicing Prediction with Discrete Compositional Energy Network"](https://dl.acm.org/doi/10.1145/3450439.3451857)

Scripts are using RNA-seq data from the ARCHS4 database. The generated dataset is available [here](https://doi.org/10.21979/N9/FFN0XH).

Download the following files necessary to generate the inputs and labels for the model to the ./data folder:

`wget "https://s3.amazonaws.com/mssm-seq-matrix/human_matrix.h5"`

`wget --timestamping 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz' -O hg38.fa.gz`

`gunzip hg38.fa.gz`

Run the main bash script for inputs, labels generation:

`generate_inputsLabels_trainTest.sh`
