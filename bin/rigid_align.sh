#!/bin/bash

# Rigid-body alignment script for small models
# Large models are best run in parallel on
# multiple processors; see the README for directions
# on how to do this, and on how to tweak options

# This script assumes the appropriate bin directory is in your path

cd bin
export PATH=$PATH:$(pwd)
#export PATH=$PATH:~/research/lib/registration/tps_alignment/bin.Linux
# set up the input files and output directory
#
# NOTE: This script does not preprocess meshes; large
# models should do so for a big performance gain.
mkdir -p global
cd global
\ls ../*.ply > ply_list.txt

# run pairwise correspondences without locally-weighted ICP
# remove the -nolocal parameter to enable locally-weighted ICP
# remove -read_xf to ignore any .xf files in the directory
#
# The sample rate is appropriate for small models; it should be
# lowered for large models
echo running pairwise correspondence...
for i in ../*.ply
  do (time correspond $i -kd 5 -read_xf -nolocal ply_list.txt $1 > corr_$(basename $i ply)txt) 2> err_$(basename $i ply)txt
done

# construct the input file for global registration
echo construct input file for global registration...
nscans=$(\ls ../*.ply | wc -w)
echo $nscans $nscans > global_in.txt
\ls ../*.ply >> global_in.txt
\ls corr*.txt >> global_in.txt
mkdir -p points rigid

# run global registration, writing a .xf file for each rigid transform in global/rigid
# remove -read_xf to ignore existing .xf files
# remove -write_xf to write a transformed mesh instead of a .xf
# add -nr nonrigid/ to write non-rigid alignments (as meshes) to global/nonrigid
echo running global registration...
(time global_reg -read_xf -write_xf -min_target_dist 3 -r rigid/ < global_in.txt > global_out.txt) 2> global_err.txt

# add symlinks to the input meshes in the same directory as the output .xf files
# NOTE: comment this out if you remove -write_xf!
#echo creating transformations...
#cd rigid
#for i in *.xf; do ln -s ../../${i%xf}ply ${i%xf}ply; done

#uncomment if you want to vrip
# rm -rf vrip
# mkdir -p vrip
# cd vrip
# for i in ../*.xf; do mesh_filter ../../../$(basename $i xf)ply -xform $i -erode -erode -erode nogrid:$(basename $i xf)ply; done
# for i in *.ply; do echo "bmesh $i 0 0 0 0 0 0 1" >> vrip.conf; done
# autorient vrip.conf
# time vripnew vrip.vri vrip.conf vrip.conf 0.25 -ramplength 3
# time vripsurf vrip.vri out.ply
# popd
# done
