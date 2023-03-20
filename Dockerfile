#LABEL maintainer="Eléonore Schneegans and Nurun Fancy

## Use rstudio installs binaries from RStudio's RSPM service by default,
## Uses the latest stable ubuntu, R and Bioconductor versions. Created on unbuntu 20.04, R 4.0 and BiocManager 3.12
FROM rocker/rstudio:4.2.2
#FROM r-base:4.0.2

## Add packages dependencies
RUN apt-get update \
	&& apt-get install -y --no-install-recommends apt-utils \
	&& apt-get install -y --no-install-recommends \
	## Basic deps
	gdb \
	libxml2-dev \
	python3-pip \
	python3 \
  python3-setuptools \
  python3-dev \
	libz-dev \
	liblzma-dev \
	libbz2-dev \
	libpng-dev \
	libgit2-dev \
	## sys deps from bioc_full
	pkg-config \
	fortran77-compiler \
	byacc \
	automake \
	curl \
	## This section installs libraries
	libpcre2-dev \
	libnetcdf-dev \
	libhdf5-serial-dev \
	libfftw3-dev \
	libopenbabel-dev \
	libopenmpi-dev \
	libxt-dev \
	libudunits2-dev \
	libgeos-dev \
	libproj-dev \
	libcairo2-dev \
	libtiff5-dev \
	libreadline-dev \
	libgsl0-dev \
	libgslcblas0 \
	libgtk2.0-dev \
	libgl1-mesa-dev \
	libglu1-mesa-dev \
	libgmp3-dev \
	libhdf5-dev \
	libncurses-dev \
	libbz2-dev \
	libxpm-dev \
	liblapack-dev \
	libv8-dev \
	libgtkmm-2.4-dev \
	libmpfr-dev \
	libmodule-build-perl \
	libapparmor-dev \
	libprotoc-dev \
	librdf0-dev \
	libmagick++-dev \
	libsasl2-dev \
	libpoppler-cpp-dev \
	libprotobuf-dev \
	libpq-dev \
	libperl-dev \
	## software - perl extentions and modules
	libarchive-extract-perl \
	libfile-copy-recursive-perl \
	libcgi-pm-perl \
	libdbi-perl \
	libdbd-mysql-perl \
	libxml-simple-perl \
	libmysqlclient-dev \
	default-libmysqlclient-dev \
	libgdal-dev \
	## new libs
	libglpk-dev \
	## Databases and other software
	sqlite \
	openmpi-bin \
	mpi-default-bin \
	openmpi-common \
	openmpi-doc \
	tcl8.6-dev \
	tk-dev \
	default-jdk \
	#imagemagick \
	#tabix \
	#ggobi \
	#graphviz \
	#protobuf-compiler \
	#jags \
	## Additional resources
	xfonts-100dpi \
	xfonts-75dpi \
	biber \
	libsbml5-dev \
	## qpdf needed to stop R CMD Check warning
	qpdf \
	&& apt-get clean \
	&& rm -rf /var/lib/apt/lists/*

RUN apt-get update \
    && apt-get install -y libcurl4-openssl-dev
RUN apt-get update \
    && apt-get install -y libcairo2-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev libxt-dev libharfbuzz-dev libfribidi-dev

# Install mofapy2
RUN python3 -m pip install 'https://github.com/bioFAM/mofapy2/tarball/master'

#RUN echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] http://packages.cloud.google.com/apt cloud-sdk main" | \
#tee -a /etc/apt/sources.list.d/google-cloud-sdk.list \
#&& curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | \
#tee /usr/share/keyrings/cloud.google.gpg && apt-get update -y \
#&& apt-get install google-cloud-sdk -y

#RUN curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" \
#-o "awscliv2.zip"

#RUN unzip awscliv2.zip && ./aws/install \
#&& rm -rf awscliv2.zip

RUN install2.r -e \
      cli \
      statmod \
      assertthat \
      BiocManager \
      devtools \
      remotes \
      enrichR \
      magrittr \
      stats \
      lme4 \
      dplyr \
      matrixStats \
      purrr \
      rmarkdown \
      tidyverse \
      paletteer \
      data.table \
      ggpubr \
      igraph \
      ggplot2 \
      cowplot \
      httr \
      jsonlite \
      reshape2 \
      ghql \
      viridis \
      tidyr \
      tibble \
      visNetwork \
      ggpubr \
      viridis \
      ggpubr \
      RColorBrewer \
      corrplot \
      stringr \
      purrr \
      ggrepel \
      SNFtool \
      && rm -rf /tmp/downloaded_packages

## Install Bioconductor packages
COPY ./misc/requirements-bioc.R .
RUN apt-get update \
 && apt-get install -y --no-install-recommends \
   libfftw3-dev \
   gcc && apt-get clean \
 && rm -rf /var/lib/apt/lists/*
RUN Rscript -e 'requireNamespace("BiocManager"); BiocManager::install(ask=F);' \
&& Rscript requirements-bioc.R \
&& rm -rf /tmp/downloaded_packages

RUN install2.r -e \
      ClassDiscovery \
      && rm -rf /tmp/downloaded_packages

# Install bioconductor dependencies
RUN R --vanilla -e "\
  if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager', repos = 'https://cran.r-project.org'); \
  sapply(c('rhdf5', 'dplyr', 'tidyr', 'reshape2', 'pheatmap', 'corrplot', \
           'ggplot2', 'ggbeeswarm', 'scales', 'GGally', 'doParallel', 'RColorBrewer', \
           'cowplot', 'ggrepel', 'foreach', 'reticulate', 'HDF5Array', 'DelayedArray', \
           'ggpubr', 'forcats', 'Rtsne', 'uwot', \
           'systemfonts', 'ragg', 'Cairo', 'ggrastr', 'basilisk', 'mvtnorm'), \
         BiocManager::install)"

RUN Rscript -e  'reticulate::py_config()'

## Install from GH the following
RUN installGithub.r cran/heatmap.plus \
    xlucpu/MOVICS \
&& rm -rf /tmp/downloaded_packages

## Install Omix package
# Copy description
WORKDIR Omix
ADD . .

# Run R CMD check - will fail with any errors or warnings
RUN Rscript -e "devtools::check()"
# Install R package from source
RUN Rscript -e "remotes::install_local()"
RUN rm -rf *
