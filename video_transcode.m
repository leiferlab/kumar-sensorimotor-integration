function success = video_transcode(input_video_file)
    output_video_file = [input_video_file(1:end-4), '.mkv'];
    if strcmp(computer, 'PCWIN64')
        %windows
        %ffmpeg -r 30 -i TestCameraFrames\Frame_%06d.tif -c:v hevc_nvenc -rc:v vbr_hq -preset:v slow CameraFrames.mp4
        system([fileparts(mfilename('fullpath')), filesep, 'LabviewVIs', filesep, 'ffmpeg -y -i ', input_video_file, ' -c:v libxvid -q:v 1 ', output_video_file])
    else
        %linux
        system(['ffmpeg -y -i ', input_video_file, ' -c:v libxvid -q:v 1 ', output_video_file]);
    end
    if exist(input_video_file, 'file') == 2 && exist(output_video_file, 'file') == 2
        delete(input_video_file);
        success = true;
    else
        success = false;
    end
end