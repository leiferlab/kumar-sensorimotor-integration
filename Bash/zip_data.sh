#!/bin/bash
folder_name=${1%/}

#read the last line of the log, get job IX
log_name=$folder_name/status.csv
last_line=$(tac $log_name |egrep -m 1 .)
IFS=',' read -ra ADDR <<< "$last_line"
job_ID=${ADDR[4]}

UpdateLog.sh $folder_name zip_data $job_ID START Zipping_Data_And_Deleting_Raw_Images


#get how many img files are there in the CameraFrames folder
#camera_img_count=$(ls $folder_name/CameraFrames/*.png | wc -l)
#if [ "$camera_img_count" -gt "0" ]; then 
#	# img file exists, zip them
#	cd $folder_name
#	zip -m -u raw_images.zip CameraFrames/*.png > /dev/null
#fi

#delete uncompressed camera files
if test -f $folder_name/CameraFrames.mkv; then 
	camera_img_count=$(ls $folder_name/DecodedCameraFrames/*.png | wc -l)
	rm -rf $folder_name/DecodedCameraFrames
fi

#get how many img files are there in the ConvertedProjectorFrames folder
projector_img_count=$(ls $folder_name/ConvertedProjectorFrames/*.png | wc -l)
if [ "$projector_img_count" -gt "0" ]; then 
	# img file exists, zip them
	cd $folder_name
	zip -m -u raw_images.zip ConvertedProjectorFrames/*.png > /dev/null
fi
UpdateLog.sh $folder_name zip_data $job_ID COMPLETE $camera_img_count'_camera_files_deleted_'$projector_img_count'_projector_img_files_zipped'
