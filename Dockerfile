#LABEL maintainer="Eléonore Schneegans and Nurun Fancy

## Uses the latest stable ubuntu, R and Bioconductor versions.
## Created on unbuntu 20.04, R 4.0 and BiocManager 3.12

FROM rocker/rstudio:4.2.3

RUN apt-get update \
	&& apt-get install -y --no-install-recommends apt-utils \
	&& apt-get install -y --no-install-recommends \
	python3 \
	python3-setuptools \
	python3-dev \
	python3-pip \
	libcurl4-openssl-dev \
	libcairo2-dev \
	libfreetype6-dev \
	libpng-dev \
	libtiff5-dev \
	libjpeg-dev \
	libxt-dev \
	libharfbuzz-dev \
	libfribidi-dev \
	## Basic deps
	gcc \
	gdb \
	libxml2-dev \
	libz-dev \
	liblzma-dev \
	libbz2-dev \
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
	libudunits2-dev \
	libgeos-dev \
	libproj-dev \
	libreadline-dev \
	libgsl0-dev \
	libgslcblas0 \
	libgtk2.0-dev \
	libgl1-mesa-dev \
	libglu1-mesa-dev \
	libgmp3-dev \
	libhdf5-dev \
	libncurses-dev \
	libxpm-dev \
	liblapack-dev \
	libv8-dev \
	libgtkmm-2.4-dev \
	libmpfr-dev \
	libmodule-build-perl \
	libapparmor-dev \
	libprotoc-dev \
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
	#cmake
	build-essential \
	wget \
	&& apt-get clean \
	&& rm -rf /var/lib/apt/lists/* \
	&& wget https://github.com/Kitware/CMake/releases/download/v3.24.1/cmake-3.24.1-Linux-x86_64.sh \
      -q -O /tmp/cmake-install.sh \
      && chmod u+x /tmp/cmake-install.sh \
      && mkdir /opt/cmake-3.24.1 \
      && /tmp/cmake-install.sh --skip-license --prefix=/opt/cmake-3.24.1 \
      && rm /tmp/cmake-install.sh \
      && ln -s /opt/cmake-3.24.1/bin/* /usr/local/bin \
  # Install mofapy2
  && python3 -m pip install mofapy2

#Set CRAN mirror
RUN echo 'options(repos = c(CRAN = "https://cloud.r-project.org"))' \
>>"${R_HOME}/etc/Rprofile.site" \

#Install CRAN pkgs
&& install2.r -e \
BiocManager \
basetheme \
data.table \
devtools \
doParallel \
DT \
dplyr \
enrichR \
foreach \
forcats \
ggplot2 \
ggpubr \
ggrastr \
ggrepel \
ggsignif \
GGally \
ggbeeswarm \
ghql \
httr \
igraph \
jsonlite \
lme4 \
magrittr \
Matrix \
matrixStats \
mvtnorm \
paletteer \
pheatmap \
plotly \
purrr \
psych \
ragg \
rlang \
RColorBrewer \
readr \
remotes \
reshape2 \
reticulate \
rmarkdown \
tibble \
tidyr \
stringr \
scales \
SNFtool \
statmod \
systemfonts \
viridis \
visNetwork \
IntNMF \
ActivePathways \
&& rm -rf /tmp/downloaded_packages

## Install Bioconductor packages

COPY ./misc/requirements-bioc.R .

RUN Rscript -e 'requireNamespace("BiocManager"); BiocManager::install(ask=F);' \
&& Rscript requirements-bioc.R \
&& rm -rf /tmp/downloaded_packages \
&& install2.r -e ClassDiscovery

# RUN Rscript -e  'reticulate::py_config()'

WORKDIR Omix
ADD . .

# Run R CMD check - will fail with any errors or warnings
RUN Rscript -e "devtools::check()"

# Install R package from source
RUN Rscript -e "remotes::install_local()"
RUN rm -rf *
