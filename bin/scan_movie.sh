#!/bin/bash

# format: ./scan_movie base_file

base_file=$1

for f in `ls ${base_file}*.png`
do
	f_jpg=${f%png}jpg
	convert ${f} -quality 100 ${f_jpg}
done

rm ${base_file}*.png

ffmpeg -f image2 -r 20 -b 800000 -i ${base_file}%05d.jpg ${base_file}.avi

rm ${base_file}*.jpg
