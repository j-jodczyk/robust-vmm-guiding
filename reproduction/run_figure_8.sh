#!/bin/bash

# run this script from within the mitsuba directory

scene_dir=scenes

# set LD_LIBRARY_PATH
export LD_LIBRARY_PATH=./dist:${LD_LIBRARY_PATH}

# run experiments
echo scene\;spp_pt\;relMSE_pt\;spp_Li\;relMSE_Li\;spp_parallax\;relMSE_parallax\;spp_adaptive\;relMSE_adaptive\;spp_cos\;relMSE_cos\;spp_BSDF\;relMSE_BSDF\;spp_NEE\;relMSE_NEE

echo Running PT baseline in living_room scene ... 1>&2
./dist/mitsuba -z ${scene_dir}/living_room/pg_script.xml -D parallaxCompensation=false -D splitAndMerge=false -D useCosineProduct=false -D useBSDFProduct=false -D bsdfProbability=0.5  -D useNee=false -D trainingSamples=0   -D sampleCount=2960 -o ${scene_dir}/living_room/living_room_pt.exr > ${scene_dir}/living_room/living_room_pt.log
echo Running Li in living_room scene ... 1>&2
./dist/mitsuba -z ${scene_dir}/living_room/pg_script.xml -D parallaxCompensation=false -D splitAndMerge=false -D useCosineProduct=false -D useBSDFProduct=false -D bsdfProbability=0.5  -D useNee=false -D trainingSamples=832 -D sampleCount=2336 -o ${scene_dir}/living_room/living_room_pg_li.exr > ${scene_dir}/living_room/living_room_pg_li.log
echo Running +parallax in living_room scene ... 1>&2
./dist/mitsuba -z ${scene_dir}/living_room/pg_script.xml -D parallaxCompensation=true  -D splitAndMerge=false -D useCosineProduct=false -D useBSDFProduct=false -D bsdfProbability=0.5  -D useNee=false -D trainingSamples=668 -D sampleCount=2016 -o ${scene_dir}/living_room/living_room_pg_li_parallax.exr > ${scene_dir}/living_room/living_room_pg_li_parallax.log
echo Running +adaptive in living_room scene ... 1>&2
./dist/mitsuba -z ${scene_dir}/living_room/pg_script.xml -D parallaxCompensation=true  -D splitAndMerge=true  -D useCosineProduct=false -D useBSDFProduct=false -D bsdfProbability=0.5  -D useNee=false -D trainingSamples=588 -D sampleCount=1984 -o ${scene_dir}/living_room/living_room_pg_li_parallax_ssmem.exr > ${scene_dir}/living_room/living_room_pg_li_parallax_ssmem.log
echo Running +cos in living_room scene ... 1>&2
./dist/mitsuba -z ${scene_dir}/living_room/pg_script.xml -D parallaxCompensation=true  -D splitAndMerge=true  -D useCosineProduct=true  -D useBSDFProduct=false -D bsdfProbability=0.5  -D useNee=false -D trainingSamples=428 -D sampleCount=1472 -o ${scene_dir}/living_room/living_room_pg_li_parallax_ssmem_cos.exr > ${scene_dir}/living_room/living_room_pg_li_parallax_ssmem_cos.log
echo Running +BSDF in living_room scene ... 1>&2
./dist/mitsuba -z ${scene_dir}/living_room/pg_script.xml -D parallaxCompensation=true  -D splitAndMerge=true  -D useCosineProduct=true  -D useBSDFProduct=true  -D bsdfProbability=0.25 -D useNee=false -D trainingSamples=384 -D sampleCount=1440 -o ${scene_dir}/living_room/living_room_pg_li_parallax_ssmem_cos_bsdf.exr > ${scene_dir}/living_room/living_room_pg_li_parallax_ssmem_cos_bsdf.log
echo Running +NEE in living_room scene ... 1>&2
./dist/mitsuba -z ${scene_dir}/living_room/pg_script.xml -D parallaxCompensation=true  -D splitAndMerge=true  -D useCosineProduct=true  -D useBSDFProduct=true  -D bsdfProbability=0.25 -D useNee=true  -D trainingSamples=232 -D sampleCount=928  -o ${scene_dir}/living_room/living_room_pg_li_parallax_ssmem_cos_bsdf_nee.exr > ${scene_dir}/living_room/living_room_pg_li_parallax_ssmem_cos_bsdf_nee.log

echo -n living_room\;
python3 reproduction/stats_figure_8.py ${scene_dir}/living_room/living_room_pt.exr ${scene_dir}/living_room/living_room_pt.log \
                                       ${scene_dir}/living_room/living_room_pg_li.exr ${scene_dir}/living_room/living_room_pg_li.log \
                                       ${scene_dir}/living_room/living_room_pg_li_parallax.exr ${scene_dir}/living_room/living_room_pg_li_parallax.log \
                                       ${scene_dir}/living_room/living_room_pg_li_parallax_ssmem.exr ${scene_dir}/living_room/living_room_pg_li_parallax_ssmem.log \
                                       ${scene_dir}/living_room/living_room_pg_li_parallax_ssmem_cos.exr ${scene_dir}/living_room/living_room_pg_li_parallax_ssmem_cos.log \
                                       ${scene_dir}/living_room/living_room_pg_li_parallax_ssmem_cos_bsdf.exr ${scene_dir}/living_room/living_room_pg_li_parallax_ssmem_cos_bsdf.log \
                                       ${scene_dir}/living_room/living_room_pg_li_parallax_ssmem_cos_bsdf_nee.exr ${scene_dir}/living_room/living_room_pg_li_parallax_ssmem_cos_bsdf_nee.log \
                                       ${scene_dir}/living_room/living_room_ref.exr

echo Running PT baseline in jewelry scene ... 1>&2
./dist/mitsuba -z ${scene_dir}/jewelry/pg_script.xml -D parallaxCompensation=false -D splitAndMerge=false -D useCosineProduct=false -D useBSDFProduct=false -D bsdfProbability=0.5  -D useNee=false -D trainingSamples=0    -D sampleCount=5148 -o ${scene_dir}/jewelry/jewelry_pt.exr > ${scene_dir}/jewelry/jewelry_pt.log
echo Running Li in jewelry scene ... 1>&2
./dist/mitsuba -z ${scene_dir}/jewelry/pg_script.xml -D parallaxCompensation=false -D splitAndMerge=false -D useCosineProduct=false -D useBSDFProduct=false -D bsdfProbability=0.5  -D useNee=false -D trainingSamples=1160 -D sampleCount=3648 -o ${scene_dir}/jewelry/jewelry_pg_li.exr > ${scene_dir}/jewelry/jewelry_pg_li.log
echo Running +parallax in jewelry scene ... 1>&2
./dist/mitsuba -z ${scene_dir}/jewelry/pg_script.xml -D parallaxCompensation=true  -D splitAndMerge=false -D useCosineProduct=false -D useBSDFProduct=false -D bsdfProbability=0.5  -D useNee=false -D trainingSamples=1008 -D sampleCount=3168 -o ${scene_dir}/jewelry/jewelry_pg_li_parallax.exr > ${scene_dir}/jewelry/jewelry_pg_li_parallax.log
echo Running +adaptive in jewelry scene ... 1>&2
./dist/mitsuba -z ${scene_dir}/jewelry/pg_script.xml -D parallaxCompensation=true  -D splitAndMerge=true  -D useCosineProduct=false -D useBSDFProduct=false -D bsdfProbability=0.5  -D useNee=false -D trainingSamples=888  -D sampleCount=3168 -o ${scene_dir}/jewelry/jewelry_pg_li_parallax_ssmem.exr > ${scene_dir}/jewelry/jewelry_pg_li_parallax_ssmem.log
echo Running +cos in jewelry scene ... 1>&2
./dist/mitsuba -z ${scene_dir}/jewelry/pg_script.xml -D parallaxCompensation=true  -D splitAndMerge=true  -D useCosineProduct=true  -D useBSDFProduct=false -D bsdfProbability=0.5  -D useNee=false -D trainingSamples=776  -D sampleCount=2624 -o ${scene_dir}/jewelry/jewelry_pg_li_parallax_ssmem_cos.exr > ${scene_dir}/jewelry/jewelry_pg_li_parallax_ssmem_cos.log
echo Running +BSDF in jewelry scene ... 1>&2
./dist/mitsuba -z ${scene_dir}/jewelry/pg_script.xml -D parallaxCompensation=true  -D splitAndMerge=true  -D useCosineProduct=true  -D useBSDFProduct=true  -D bsdfProbability=0.25 -D useNee=false -D trainingSamples=592  -D sampleCount=1952 -o ${scene_dir}/jewelry/jewelry_pg_li_parallax_ssmem_cos_bsdf.exr > ${scene_dir}/jewelry/jewelry_pg_li_parallax_ssmem_cos_bsdf.log
echo Running +NEE in jewelry scene ... 1>&2
./dist/mitsuba -z ${scene_dir}/jewelry/pg_script.xml -D parallaxCompensation=true  -D splitAndMerge=true  -D useCosineProduct=true  -D useBSDFProduct=true  -D bsdfProbability=0.25 -D useNee=true  -D trainingSamples=436  -D sampleCount=1536 -o ${scene_dir}/jewelry/jewelry_pg_li_parallax_ssmem_cos_bsdf_nee.exr > ${scene_dir}/jewelry/jewelry_pg_li_parallax_ssmem_cos_bsdf_nee.log

echo -n jewelry\;
python3 reproduction/stats_figure_8.py ${scene_dir}/jewelry/jewelry_pt.exr ${scene_dir}/jewelry/jewelry_pt.log \
                                      ${scene_dir}/jewelry/jewelry_pg_li.exr ${scene_dir}/jewelry/jewelry_pg_li.log \
                                      ${scene_dir}/jewelry/jewelry_pg_li_parallax.exr ${scene_dir}/jewelry/jewelry_pg_li_parallax.log \
                                      ${scene_dir}/jewelry/jewelry_pg_li_parallax_ssmem.exr ${scene_dir}/jewelry/jewelry_pg_li_parallax_ssmem.log \
                                      ${scene_dir}/jewelry/jewelry_pg_li_parallax_ssmem_cos.exr ${scene_dir}/jewelry/jewelry_pg_li_parallax_ssmem_cos.log \
                                      ${scene_dir}/jewelry/jewelry_pg_li_parallax_ssmem_cos_bsdf.exr ${scene_dir}/jewelry/jewelry_pg_li_parallax_ssmem_cos_bsdf.log \
                                      ${scene_dir}/jewelry/jewelry_pg_li_parallax_ssmem_cos_bsdf_nee.exr ${scene_dir}/jewelry/jewelry_pg_li_parallax_ssmem_cos_bsdf_nee.log \
                                      ${scene_dir}/jewelry/jewelry_ref.exr
