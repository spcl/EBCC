# Use a manylinux base image
FROM quay.io/pypa/manylinux_2_28_x86_64 as base

# Set environment variables
ENV LANG C.UTF-8
ENV LC_ALL C.UTF-8

# Install system dependencies
RUN dnf install -y \
    git \
    wget \
    && dnf clean all
RUN dnf install -y epel-release dnf-plugins-core redhat-rpm-config && \
    dnf --enablerepo=powertools install -y netcdf-devel-4.7.0-3.el8.x86_64

# Set working directory
WORKDIR /workspace

FROM base as build_plugin

# Install Python dependencies (example: using requirements.txt)
# RUN /opt/python/cp39-cp39/bin/pip install --no-cache-dir -r requirements.txt

RUN wget https://github.com/HDFGroup/hdf5/archive/refs/tags/hdf5_1.14.4.3.tar.gz && tar -xzf hdf5_1.14.4.3.tar.gz && rm -f hdf5_1.14.4.3.tar.gz
RUN mv hdf5-hdf5_1.14.4.3/ hdf5-1.14.4-3

RUN cp -r hdf5-1.14.4-3/config/cmake/scripts/* .
RUN ctest -S HDF5config.cmake,BUILD_GENERATOR=Unix -C Release -VV -O hdf5.log
RUN sed -i 's/cpack_skip_license=FALSE/cpack_skip_license=TRUE/' build/HDF5-1.14.4.3-Linux.sh
RUN build/HDF5-1.14.4.3-Linux.sh
RUN wget https://github.com/facebook/zstd/archive/refs/tags/v1.5.6.tar.gz && tar -xvf v1.5.6.tar.gz && rm -f v1.5.6.tar.gz \
    && cd zstd-1.5.6 && make
RUN wget https://github.com/uclouvain/openjpeg/archive/refs/tags/v2.5.2.tar.gz && tar -xvf v2.5.2.tar.gz && rm -f v2.5.2.tar.gz \
    && cd openjpeg-2.5.2 && mkdir build && cd build && cmake .. -DCMAKE_BUILD_TYPE=Release && make
RUN git clone https://github.com/rafat/wavelib.git && cd wavelib && git reset --hard a92456d2e20451772dd76c2a0a3368537ee94184 \
    && mkdir build && cd build && cmake -DCMAKE_POSITION_INDEPENDENT_CODE=ON .. && make

RUN git clone https://github.com/luigifusco/compression-filter.git && cd compression-filter/src \
    && CPATH=/workspace/HDF_Group/HDF5/1.14.4.3/include/:/workspace/zstd-1.5.6/lib/:/workspace/openjpeg-2.5.2/src/lib/openjp2/:/workspace/openjpeg-2.5.2/build/src/lib/openjp2/:/workspace/wavelib/header \
    LIBRARY_PATH=/workspace/wavelib/build/Bin/:/workspace/openjpeg-2.5.2/build/bin/:/workspace/zstd-1.5.6/lib/:/workspace/HDF_Group/HDF5/1.14.4.3/lib gcc -Wall -g -O2 -shared -fPIC -o libh5z_j2k.so -l:libzstd.a -l:libopenjp2.a -l:libhdf5.a h5z_j2k.c

FROM base as final
RUN git clone https://github.com/luigifusco/compression-filter.git
COPY --from=build_plugin /workspace/compression-filter/src/libh5z_j2k.so /workspace/compression-filter/src

CMD ["/bin/bash"]
