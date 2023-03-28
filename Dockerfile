#LABEL maintainer="Eléonore Schneegans and Nurun Fancy

## Use rstudio installs binaries from RStudio's RSPM service by default,
## Uses the latest stable ubuntu, R and Bioconductor versions. Created on unbuntu 20.04, R 4.0 and BiocManager 3.12

FROM rocker/rstudio:4.0.2


RUN apt-get update \
	&& apt-get install -y --no-install-recommends apt-utils \
	&& apt-get install -y --no-install-recommends \
	## Basic deps
	gdb \
	libxml2-dev \
	libz-dev \
	liblzma-dev \
	libbz2-dev \
	libgit2-dev \
	python3 \
	python3-setuptools \
	python3-dev \
	python3-pip \
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
	libraptor2-dev \
	librasqal3-dev \
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
	imagemagick \
	tabix \
	ggobi \
	graphviz \
	protobuf-compiler \
	jags \
	## Additional resources
	xfonts-100dpi \
	xfonts-75dpi \
	biber \
	libsbml5-dev \
	## qpdf needed to stop R CMD Check warning
	qpdf \
	## MOFA
	libcurl4-openssl-dev \
	libcairo2-dev \
	libfreetype6-dev \
	libpng-dev \
	libtiff5-dev \
	libjpeg-dev \
	libxt-dev \
	libharfbuzz-dev \
	libfribidi-dev \
	&& apt-get clean \
	&& rm -rf /var/lib/apt/lists/*


# Install mofapy2
RUN python3 -m pip install 'https://github.com/bioFAM/mofapy2/tarball/master'

RUN install2.r -e \
      cli \
      statmod \
      assertthat \
      BiocManager \
      devtools \
      remotes \
      enrichR \
      magrittr \
      lme4 \
      dplyr \
      matrixStats \
      purrr \
      Matrix \
      rmarkdown \
      tidyverse \
      paletteer \
      data.table \
      ggpubr \
      igraph \
      plotly \
      ggplot2 \
	  ggbeeswarm \
	  GGally \
	  ggrastr \
      cowplot \
      httr \
      jsonlite \
      reshape2 \
	  pheatmap \
      ghql \
      viridis \
      tidyr \
      tibble \
      visNetwork \
      ggsignif \
      DT \
      RColorBrewer \
	  reticulate \
      corrplot \
	  doParallel \
      stringr \
      ggrepel \
	  foreach \
	  forcats \
	  systemfonts \
	  ragg \
      SNFtool \
	  Cairo \
	  mvtnorm \
	  scales \
	  ClassDiscovery\
      && rm -rf /tmp/downloaded_packages


## Install Bioconductor packages
COPY ./misc/requirements-bioc.R .

RUN Rscript -e 'requireNamespace("BiocManager"); BiocManager::install(ask=F);' \
&& Rscript requirements-bioc.R \
&& rm -rf /tmp/downloaded_packages

RUN Rscript -e  'reticulate::py_config()'

WORKDIR Omix
ADD . .

# Run R CMD check - will fail with any errors or warnings
RUN Rscript -e "devtools::check()"
# Install R package from source
RUN Rscript -e "remotes::install_local()"
RUN rm -rf *
