#
# based on
# https://github.com/benjamin-heasly/mitsuba-docker/blob/master/rgb/Dockerfile
#

# builder
FROM ubuntu:22.04 AS builder
LABEL authors="Alexander Rath <alexander.rath@dfki.de>"

ENV TZ=Europe
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

RUN apt-get update && \
    apt-get install -y \
        build-essential \
        scons \
        libgl1-mesa-dev \
        libglu1-mesa-dev \
        libxxf86vm-dev \
        libpng-dev \
        libjpeg-dev \
        libilmbase-dev \
        libxerces-c-dev \
        libboost-all-dev \
        libopenexr-dev \
        libglew-dev \
        libpcrecpp0v5 \
        libeigen3-dev \
        libfftw3-dev \
        wget && \
    apt-get clean && \
    apt-get autoclean && \
    apt-get autoremove

WORKDIR /mitsuba
COPY mitsuba .

RUN cp build/config-linux-gcc.py config.py && scons -j8

# mitsuba
FROM ubuntu:22.04 AS mitsuba
LABEL authors="Alexander Rath <alexander.rath@dfki.de>"

ENV TZ=Europe
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

RUN apt-get update && \
    apt-get install -y \
        libgl1-mesa-dev \
        libglu1-mesa-dev \
        libxxf86vm-dev \
        libpng-dev \
        libjpeg-dev \
        libilmbase-dev \
        libxerces-c-dev \
        libboost-all-dev \
        libopenexr-dev \
        libglew-dev \
        libpcrecpp0v5 \
        libeigen3-dev \
        libfftw3-dev && \
    apt-get clean && \
    apt-get autoclean && \
    apt-get autoremove

WORKDIR /mitsuba
ENV MITSUBA_DIR /mitsuba
ENV PYTHONPATH /mitsuba/dist/python:/mitsuba/dist/python/2.7:
ENV PATH /mitsuba/wrapper:/mitsuba/dist:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
ENV LD_LIBRARY_PATH /mitsuba/dist:/usr/local/lib:

COPY --from=builder /usr/local /usr/local
COPY --from=builder /mitsuba .

# Set default scene file
ENV SCENE_FILE scene/dining-room/dining-room.xml

# Use the environment variable for the scene file
CMD ["sh", "-c", "mitsuba $SCENE_FILE"]