FROM ubuntu:22.04

LABEL name="Evan Petrimoulx"
LABEL maintainer="Evan Petrimoulx"

RUN apt-get update && apt-get install -y \
    build-essential \
    gfortran \
    cmake \
    bash

COPY ./ /matrix-elements-dq
WORKDIR /matrix-elements-dq
RUN make all
