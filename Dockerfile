FROM continuumio/miniconda3:4.8.2
LABEL authors="chris.cheshire@crick.ac.uk" \
      description="Docker image containing all requirements for tobias https://github.com/loosolab/TOBIAS"

RUN apt-get update \
 && apt-get install -y --no-install-recommends \
 procps=2:3.3.15-2 \
 libz-dev \
 build-essential \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*

 RUN pip install tobias