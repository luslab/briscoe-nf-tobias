#!/bin/bash

VERSION="0.12.10"
ECHO "Building v$VERSION"

docker build -f Dockerfile -t luslab/nf-modules-tobias:$VERSION .