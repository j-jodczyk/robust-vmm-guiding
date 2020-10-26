#!/bin/bash

# run this script from within the mitsuba directory

scene_dir=scenes

run_scene() {
    echo Running experiment for the $1 scene ... 1>&2
    ./dist/mitsuba -z ${scene_dir}/$1/pg_script.xml -D useNee=false -D trainingSamples=$2 -D sampleCount=$3 -o ${scene_dir}/$1/$1_table3.exr > ${scene_dir}/$1/$1_table3.log
    echo -n $1\;
    python3 reproduction/stats_table_3.py ${scene_dir}/$1/$1_table3.exr ${scene_dir}/$1/$1_table3.log ${scene_dir}/$1/$1_ref.exr
}


# set LD_LIBRARY_PATH
export LD_LIBRARY_PATH=./dist:${LD_LIBRARY_PATH}

# run experiments
echo -n scene\;
python3 reproduction/stats_table_3.py --header

run_scene bathroom 140 316
run_scene clocks 244 4460
run_scene country_kitchen_day 208 2468
run_scene country_kitchen_night 528 2032
run_scene glossy_cbox 212 4296
run_scene jewelry 272 3496
run_scene kitchen 228 1908
run_scene kitchenette 248 708
run_scene living_room 248 1732
run_scene pool 604 2684
run_scene torus 620 4244
