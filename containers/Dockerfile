#==============================================================================
# CGS GAMERA (Serial) Containerization
# 
# This Dockerfile allows one to create a build environment for the GAMERA model.
# The base image/layer is Intel's OneAPI HPC Toolkit allowing for Intel-based
# compilers.
#
# Run:
#   docker build -t <image_tag> . && docker run --rm -it <image_tag>
#
# Note: The initial build will take longer as Docker will be retrieving the
#       base image (~5.5GB). This time will vary depending upon your network.
# Note: The base image may fail to build (known issue). If it does, that
#       means that Intel's container repository is having network issues.
#       A work-around is in development.
# 
# 2024, Brent Smith
#==============================================================================
FROM intel/oneapi-hpckit:2023.2.1-devel-ubuntu22.04
# Intel OneAPI Docker Images:
#     https://github.com/intel/oneapi-containers/tree/master/images
#     https://hub.docker.com/r/intel/oneapi-hpckit

# Directory Structure
RUN mkdir -p /cdf /hdf5 /miniconda3 /kaiju

#------------------------------------------------------------------------------
# NASA CDF Library (version 3.9.0)
ENV CDF_FILE=cdf39_0-dist-all.tar.gz
ADD https://spdf.gsfc.nasa.gov/pub/software/cdf/dist/cdf39_0/linux/${CDF_FILE} /cdf/${CDF_FILE}
RUN cd /cdf && tar xf ${CDF_FILE} && cd cdf39_0-dist &&\
    make OS=linux ENV=gnu CURSES=no all &&\
    make INSTALLDIR=/cdf install &&\
    rm ../${CDF_FILE} &&\
    . /cdf/bin/definitions.B

# HDF5 Library
ENV HDF_FILE=hdf5-1.14.1-2
ADD https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.14/hdf5-1.14.1/src/${HDF_FILE}.tar.gz /hdf5/${HDF_FILE}.tar.gz
RUN cd /hdf5 && tar zxf ${HDF_FILE}.tar.gz && cd ${HDF_FILE} &&\
    ./configure FC=ifort --prefix=/hdf5 --enable-fortran &&\
    make && make install &&\
    rm ../${HDF_FILE}.tar.gz
ENV PATH="/hdf5/bin:${PATH}" \
    HDF5_DIR=/hdf5 \
    HDF5_INCLUDE_DIRS=/hdf5/include \
    HDF5_LIBRARIES=/hdf5/lib \
    CPATH=${HDF5_INCLUDE_DIRS} \
    INCLUDE="-I${HDF5_INCLUDE_DIRS}" \
    LD_LIBRARY_PATH="/hdf5:${LD_LIBRARY_PATH}"

# Update System
# known GPG key issue: https://github.com/intel/oneapi-containers/issues/53
RUN curl -fsSL https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS-2023.PUB | gpg --dearmor > /usr/share/keyrings/intel-oneapi-archive-keyring.gpg
RUN apt -y update && \
    apt install --no-install-recommends -y \
        build-essential \
        cmake \
        git \
        git-lfs \
        libgeos-dev \
        pkg-config && \
        # vim \
    rm -rf /var/lib/apt/lists/*

# Python Environment
ENV PY_FILE=Miniconda3-latest-Linux-x86_64.sh
COPY requirements.txt /miniconda3/
ADD https://repo.anaconda.com/miniconda/${PY_FILE} /miniconda3/${PY_FILE}
RUN cd /miniconda3 && bash ${PY_FILE} -b -u -p /miniconda3 && \
    rm ${PY_FILE} && /miniconda3/bin/conda init bash && \
    . ~/.bashrc && \
    conda config --set ssl_verify false && \
    conda config --set auto_activate_base false && \
    conda update -y --all && \
    conda create -n kaiju-3.8 python=3.8 -y && \
    conda activate kaiju-3.8 && \
    conda update -y --all && \
    pip3 install --default-timeout=100 -r /miniconda3/requirements.txt

#==============================================================================

# CGS GAMERA Setup/Compilation
ENV KAIJUHOME=/kaiju \
    FC=/hdf5/bin/h5fc

# Build pre-cloned repository code on user-defined branch (development)
COPY . /kaiju/
RUN cd /kaiju/ && git remote rm origin && \
    mkdir -p ${KAIJUHOME}/build && cd ${KAIJUHOME}/build && \
    cmake -DALLOW_INVALID_COMPILERS=ON .. && \
    make all

# Setup Running ENV and entry directory
WORKDIR /kaiju/build/
ENV PATH=${KAIJUHOME}/quickstart/loop2d:${PATH}
RUN echo "conda activate kaiju-3.8" >> ~/.bashrc && \
    echo ". /cdf/bin/definitions.B" >> ~/.bashrc && \
    echo ". \${KAIJUHOME}/scripts/setupEnvironment.sh" >> ~/.bashrc && \
    echo "PATH=\${KAIJUHOME}/build/bin:\${PATH}" >> ~/.bashrc