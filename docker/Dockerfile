# RLCtools: an R library of helper functions  
# Copyright (c) 2023-Present, Ryan L. Collins and the Dana-Farber Cancer Institute.  
# Distributed under terms of the GNU GPL v2.0 License

FROM condaforge/miniforge3:24.9.2-0
MAINTAINER "Ryan Collins <Ryan_Collins@dfci.harvard.edu>"

# Aliases for convenience
RUN echo 'alias l="ls -ltrha"' >> ~/.bashrc && \
    echo 'alias less="zless"' >> ~/.bashrc

# Install git
RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy install --no-upgrade --no-install-recommends git

# Install R dependencies with conda
RUN mamba install -n base --update-all \
    -c conda-forge -c r -c bioconda \
    r-base r-bedr r-beeswarm r-desctools r-eql r-survival r-vioplot r-viridis

# Clone RLCtools repo and checkout specified hash
ARG REPO_HASH=main
RUN git clone \
    https://github.com/RCollins13/RLCtools.git \
    opt/RLCtools \
    --branch=$REPO_HASH

# Install RLCtools from source
ARG PKG_VERSION=0.1
RUN Rscript -e "install.packages(\"/opt/RLCtools/RLCtools_$PKG_VERSION.tar.gz\", type=\"source\", repos=NULL)"

# Launch bash at runtime
CMD ["/bin/bash"]
