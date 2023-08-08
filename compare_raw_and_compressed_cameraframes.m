
    %ffmpeg -i CameraFrames\Frame_%06d.tif -threads 8 -an -vcodec ffv1 -coder 1 -context 1 -g 1 -level 3 -slices 24 -slicecrc 1 -pass 1 -passlogfile my_passlog CameraFrames.mkv
    %ffmpeg -i CameraFrames\Frame_%06d.tif -threads 8 -an -vcodec ffv1 -coder 1 -context 1 -g 1 -level 3 -slices 24 -slicecrc 1 -pass 2 -passlogfile my_passlog CameraFrames.mkv
    %ffmpeg -r 30 -i TestCameraFrames\Frame_%06d.tif -c:v hevc_nvenc -rc:v vbr_hq -preset:v slow CameraFrames.mp4
    % ffmpeg -y -r 30 -i CameraFrames\Frame_%06d.tif -c:v hevc_nvenc -rc:v vbr_hq -preset:v slow CameraFrames.mkv
    %ffmpeg -i CameraFrames.mkv -compression_algo raw -pix_fmt gray -start_number 0 DecodedCameraFrames\Frame_%06d.tif
    %ffmpeg -i CameraFrames.mp4 -compression_algo raw -pix_fmt gray -start_number 0 DecodedCameraFrames\Frame_%06d.tif


%     folder_name = 'C:\Data\20200217_HeadTailContinuousGWN_test\Data20200217_115033';
%     folder_name = 'Z:\Mochi\APIData\20200217_HeadTailContinuousGWN_test\Data20200217_115033';
%     folder_name = 'Z:\Mochi\APIData\20200224_launched_from_main_test\Data20200224_155118';
    folder_name = 'C:\Data\20200224_launched_from_main_test\Data20200224_155118';
    parameters = load_parameters(folder_name); %load experiment parameters
    
    %% STEP 1: initialize %%
    number_of_images_for_median_projection = parameters.PostAnalysisNumberofMedianFilterImages;
    ImageSize = 2*round(parameters.TrackingWormBoxHalfSize * parameters.CameraPixeltommConversion);
    image_size = [ImageSize, ImageSize];
    

    %% STEP 2: Load images and labview tracks from the directory %%
    camera_image_directory = [folder_name, filesep, 'CameraFrames', filesep];
    image_files = dir([camera_image_directory, '*.tif']); %get all the image files
    
    decoded_camera_image_directory = [folder_name, filesep, 'DecodedCameraFrames', filesep];
    decoded_image_files = dir([camera_image_directory, '*.tif']); %get all the image files

    
    %% STEP 3: Get the median projection for use as background%%
    medianProj = imread([camera_image_directory, image_files(1).name]);
    medianProjCount = min(number_of_images_for_median_projection, length(image_files) - 1); 
    medianProj = zeros(size(medianProj,1), size(medianProj,2), medianProjCount);
    for frame_index = 1:medianProjCount
        image_file_name = image_files(floor((length(image_files)-1)*frame_index/medianProjCount)).name
        curImage = imread([camera_image_directory, image_file_name]);
        medianProj(:,:,frame_index) = curImage;
        imagesc(curImage)
        colorbar
        pause
        
        decoded_file_name = decoded_image_files(floor((length(decoded_image_files)-1)*frame_index/medianProjCount)).name
        curDecodedImage = imread([decoded_camera_image_directory, decoded_file_name]);
        imagesc(curDecodedImage)
        colorbar
        pause
        
        diff_image = double(curImage) - double(curDecodedImage);
        sum(diff_image(:))
        imagesc(diff_image);
        
        colorbar
        pause;
    end
    medianProj = median(medianProj, 3);
    medianProj = uint8(medianProj);
