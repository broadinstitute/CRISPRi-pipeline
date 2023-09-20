FROM us.gcr.io/broad-dsde-methods/terra-jupyter-minimal-gpu-base:0.0.2
# Build off of image https://github.com/sjfleming/terra-docker/blob/sf_minimal_base/terra-jupyter-minimal-gpu-base/CHANGELOG.md
# terra-jupyter-minimal-gpu-base image


WORKDIR /app
ENV PATH=$PATH:/app

COPY requirements.txt .
RUN pip3 install --break-system-packages -r requirements.txt

COPY CRISPRi_pipeline ./CRISPRi_pipeline

# Copy code for installing R and R packages from terra-jupyter-r image:
# https://github.com/sjfleming/terra-docker/blob/sf_minimal_base/terra-jupyter-r/README.md

# https://cran.r-project.org/bin/linux/ubuntu/README.html
RUN apt-get update \
    && apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 \
    && add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran40/' \
    && apt-get install -yq --no-install-recommends apt-transport-https \
    && apt update \
    && apt install -yq --no-install-recommends \
	apt-utils \
	libssh2-1-dev \
	libssl-dev \
	libcurl4-gnutls-dev \
	libgit2-dev \
	libxml2-dev \
	libgfortran-7-dev \
	r-base-dev \
	r-base-core \
	## This section installs libraries
	libnetcdf-dev \
	libhdf5-serial-dev \
	libfftw3-dev \
	libopenbabel-dev \
	libopenmpi-dev \
	libexempi3 \
	libgdal-dev \
	libcairo2-dev \
	libtiff5-dev \
	libgsl0-dev \
	libgtk2.0-dev \
	libgl1-mesa-dev \
	libglu1-mesa-dev \
	libgmp3-dev \
	libhdf5-dev \
	libncurses-dev \
	libxpm-dev \
	libv8-3.14-dev \
	libgtkmm-2.4-dev \
	libmpfr-dev \
	libudunits2-dev \
	libmodule-build-perl \
	libapparmor-dev \
	libgeos-dev \
	librdf0-dev \
	libmagick++-dev \
	libsasl2-dev \
	libpoppler-cpp-dev \
	libpq-dev \
	libperl-dev \
	libgfortran5 \
	## software - perl extentions and modules
	libarchive-extract-perl \
	libfile-copy-recursive-perl \
	libcgi-pm-perl \
	libdbi-perl \
	libdbd-mysql-perl \
	libxml-simple-perl \
	## Databases and other software
	sqlite \
	mpi-default-bin \
	openmpi-common \
	tcl8.5-dev \
	imagemagick \
	tabix \
	ggobi \
	graphviz \
	jags \
	## Additional resources
	xfonts-100dpi \
	xfonts-75dpi \
	biber \
	libzmq3-dev \
	libsbml5-dev \
	biber \
        ## for gpuMagic package
        ocl-icd-opencl-dev \
        ## for ChemmineOB package
        libeigen3-dev \
        ## for rawrr package
        mono-runtime \
    && ln -s /usr/lib/gcc/x86_64-linux-gnu/7/libgfortran.so /usr/lib/x86_64-linux-gnu/libgfortran.so \
    && ln -s /usr/lib/gcc/x86_64-linux-gnu/7/libstdc++.so /usr/lib/x86_64-linux-gnu/libstdc++.so \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*


RUN R -e 'install.packages("BiocManager")' \
    ## check version
    && R -e 'BiocManager::install(version="3.17", ask=FALSE)' \
    && R -e 'BiocManager::install(c( \
    "boot", \
    "class", \
    "cluster", \
    "codetools", \
    "foreign", \
    "kernsmooth", \
    "lattice", \
    "mass", \
    "Matrix", \
    "mgcv", \
    "nlme", \
    "nnet", \
    "rpart", \
    "Seurat", \
    "spatial", \
    "survival", \
    # GCP essentials
    "bigrquery",  \
    "googleCloudStorageR", \
    # User oriented packages
    "reticulate", \
    "remotes", \
    "devtools", \
    "tidyverse", \
    "pbdZMQ", \
    "uuid", \
    # extra packages from CRISPR-de pipeline
    "reshape2", \
    "gridExtra", \
    "grid", \
    "patchwork" \
    "SeuratObject" \
    "dplyr", \
    "tidyr"))' \
    && R -e 'BiocManager::install("DataBiosphere/Ronaldo")'

# Install Bioconductor packages found at:
# https://raw.githubusercontent.com/anvilproject/anvil-docker/master/anvil-rstudio-bioconductor/install.R
RUN R -e 'BiocManager::install(c( \
    "AnVIL", \
    "SingleCellExperiment", \
    "GenomicFeatures", \
    "GenomicAlignments", \
    "ShortRead", \
    "DESeq2", \
    "AnnotationHub", \
    "ExperimentHub", \
    "ensembldb", \
    "scRNAseq", \
    "scran", \
    "Rtsne"))'

# potential package needed for MAST?
RUN Rscript -e 'devtools::install_github(repo = "mojaveazure/loomR", ref = "develop")'
