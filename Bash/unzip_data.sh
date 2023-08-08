#!/bin/bash
folder_name=${1%/}

#read the last line of the log, get job IX
log_name=$folder_name/status.csv
last_line=$(tac $log_name |egrep -m 1 .)
IFS=',' read -ra ADDR <<< "$last_line"
job_ID=${ADDR[4]}

UpdateLog.sh $folder_name unzip_data $job_ID START Unzipping_Data

#decompress using ffmpeg
if test -f $folder_name/CameraFrames.mkv; then 
	rm -rf $folder_name/DecodedCameraFrames
	mkdir $folder_name/DecodedCameraFrames
	ffmpeg -i $folder_name/CameraFrames.mkv -pix_fmt gray -start_number 0 $folder_name/DecodedCameraFrames/Frame_%06d.png
	camera_img_count=$(ls $folder_name/DecodedCameraFrames/*.png | wc -l)
fi

if test -f $folder_name/raw_images.zip; then 
	# raw images zip exists, extract
	cd $folder_name
	proj_img_count=$(($(unzip -u raw_images.zip | wc -l) - 2))
fi

UpdateLog.sh $folder_name unzip_data $job_ID COMPLETE $camera_img_count'_camera_image_decompressed_'$proj_img_count'_projector_img_files_unzipped'
