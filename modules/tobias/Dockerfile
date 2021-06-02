FROM nfcore/base:1.14
LABEL authors="chris.cheshire@crick.ac.uk" \
      description="Docker image containing all requirements for tobias https://github.com/loosolab/TOBIAS"

RUN apt-get update \
 && apt-get install -y --no-install-recommends \
 procps=2:3.3.15-2 \
 libz-dev \
 build-essential \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*

RUN mkdir -p /.config/matplotlib && chmod 777 /.config/matplotlib
RUN mkdir -p /.cache/matplotlib && chmod 777 /.cache/matplotlib

# Update conda
RUN conda update -n base -c defaults conda

# Install conda packages
COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/nfcore-module-tobias/bin:$PATH

# Install dev fork
# WORKDIR /home
# RUN git clone https://github.com/chris-cheshire/TOBIAS.git
# WORKDIR /home/TOBIAS
# RUN /opt/conda/envs/nfcore-module-tobias/bin/pip install .
# WORKDIR /home
# RUN conda install -n nfcore-module-tobias -c conda-forge filprofiler