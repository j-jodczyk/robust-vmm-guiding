#!/bin/bash

# make sure everything is up to date
git pull
git submodule update --init --recursive

# install prerequisites (tested on Ubuntu 20.04)
# skip if there already exists a compiled version
if [ ! -f "dist/mitsuba" ]
then
    # update repository cache
    sudo apt-get update

    # core build environment
    sudo apt-get install -y build-essential scons
    # image libraries
    sudo apt-get install -y libpng-dev libjpeg-dev libilmbase-dev libopenexr-dev
    # math libraries
    sudo apt-get install -y libeigen3-dev libfftw3-dev
    # boost
    sudo apt-get install -y libboost-all-dev
    # xml library
    sudo apt-get install -y libxerces-c-dev
    # regexp library
    sudo apt-get install -y libpcrecpp0v5

    ## GUI stuff - not really needed to reproduce the results
    # graphics libraries
    sudo apt-get install -y libxxf86vm-dev libglew-dev
    # qt
    sudo apt-get install -y qt5-default qt5-qmake libqt5xmlpatterns5-dev

    ## utilities
    sudo apt-get install -y wget unzip

    ## python 3 and modules required for reproduction
    sudo apt-get install -y python3 python3-pip
    sudo -H python3 -m pip install --upgrade pip
    sudo -H python3 -m pip install numpy OpenEXR Imath
fi

# compile
scons --cfg=build/config-linux-gcc.py -j8
if [ $? -ne 0 ]
then
    echo Compilation error - please report back to the authors.
    exit 1
fi

# fetch scenes
if [ ! -d "scenes" ]
then
    # these scenes are only provided for reproduction purposes - not all scenes' licences allow for redistribution
    [ -f "scenes-reproduction.zip" ] || wget http://www.graphics-data.uni-tuebingen.de/ruppert/robust_fitting_parallax_guiding_2020/scenes-reproduction.zip
    if [ $? -eq 0 -a "$(md5sum scenes-reproduction.zip | cut -f1 -d' ')" = "f736a16efcae3df39c0af3cf3e67914d" ]
    then
        unzip scenes-reproduction.zip -d scenes
    else
        echo Failed to fetch the scene files - please report back to the authors.
        exit 1
    fi
fi

echo
echo Experiments will be run shortly.
echo Note that we cleaned up the code and made some small optimizations.
echo The code should now run slightly faster than the reported numbers when run on the same hardware.
echo Our reference hardware is using dual Intel\(R\) Xeon\(R\) Gold 5115 CPUs @ 2.40GHz \(40 threads\).
echo Depending on your hardware, runtimes will differ - we fixed the sample counts used in the experiments based on the reported numbers.
echo
echo Note that the glossy_cbox scene changed due to a fixed bug in Mitsuba. The results will differ slightly as a result of that.
echo Finally, keep in mind that the random nature of path tracing will introduce some variance to the results, especially in terms of the reported relMSE.

echo
echo Running experiments for table 3 / figure 9 \(Li cos only\) ... \(about 2 hours total / 10 minutes per scene on similar hardware\)
echo Results will be created in the scenes folder as scene_table3.exr and scene_table3.log
echo The resulting table will be stored in table3.csv \(speedup = relMSE baseline / relMSE of the result\)
sh reproduction/run_table_3.sh > table3.csv

echo
echo Running experiments for figure 8 \(about 2.5 hours total / 10 minutes per scene and configuration on similar hardware\)
echo Results will be created in the scene folders as scene_pt.exr, scene_pg_li.exr, scene_pg_li_parallax.exr,
echo scene_pg_li_parallax_ssmem_cos.exr, scene_pg_li_parallax_ssmem_cos_bsdf.exr, scene_pg_li_parallax_ssmem_cos_bsdf_nee.exr.
echo The resulting table will be stored in figure8.csv \(speedup = relMSE for Li / relMSE of the result\)
sh reproduction/run_figure_8.sh > figure8.csv
