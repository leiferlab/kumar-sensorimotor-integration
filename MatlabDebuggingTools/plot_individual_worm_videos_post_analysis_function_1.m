function plot_individual_worm_videos_post_analysis_function_1(main_folder,output_video_frame_rate,save_excel_file,path_history,demo_mode,parameters,count_number_of_reversal_trial_rails,count_number_of_reversal_trial_turns,stimulus_intensities,test_stimulus_duration,current_tracks,fps,behavior_colors,stim_color)

close all
  
current_rails_dur=test_stimulus_duration(1)*fps;

figure1=figure('Position', [40 40 820 500]);

addpath(fullfile(main_folder,'/individual_worm_imgs'));
addpath(fullfile(main_folder,'/ConvertedProjectorFrames/'));

loaded_variable=load(fullfile(main_folder,'timestamps.mat'));
processed_decoded_camera_frames = loaded_variable.processed_decoded_camera_frames;

ImageSize = 2*round(parameters.TrackingWormBoxHalfSize * parameters.CameraPixeltommConversion);
image_size = [ImageSize, ImageSize];

excel_filename= "excel_file_"  + main_folder(end-14:end) + ".csv";
video_folder_name1="individual_worm_videos";

if contains(main_folder, 'Full')
    worm_frame_matrix=count_number_of_reversal_trial_rails;
    if save_excel_file==1
       writematrix(worm_frame_matrix,fullfile(main_folder,excel_filename)) 
    end
    if isempty(stimulus_intensities)
        stimulus_intensities=0;
    end
    video_folder_name2="rails_"+num2str(ceil(stimulus_intensities(1)))+"uW_"+num2str(test_stimulus_duration(1))+"_sec";
    video_folder_name=fullfile(video_folder_name1,video_folder_name2);
elseif contains(main_folder,'Turn')
    worm_frame_matrix=count_number_of_reversal_trial_turns;
    if save_excel_file==1
       writematrix(worm_frame_matrix,fullfile(main_folder,excel_filename)) 
    end 
    if isempty(stimulus_intensities)
        stimulus_intensities=0;
    end
    video_folder_name2="turns_"+num2str(ceil(stimulus_intensities(1)))+"uW_"+num2str(test_stimulus_duration(1))+"_sec";
    video_folder_name=fullfile(video_folder_name1,video_folder_name2);
    elseif contains(main_folder,'Reversing')
    worm_frame_matrix=count_number_of_reversal_trial_turns;
    if save_excel_file==1
       writematrix(worm_frame_matrix,fullfile(main_folder,excel_filename)) 
    end
    if isempty(stimulus_intensities)
        stimulus_intensities=0;
    end
    video_folder_name2="reverse"+num2str(ceil(stimulus_intensities(1)))+"uW_"+num2str(test_stimulus_duration(1))+"_sec";
    video_folder_name=fullfile(video_folder_name1,video_folder_name2);
end

for num_iteration=1:size(worm_frame_matrix,1) 

% % %%%%%% save video only when certain criteria is met
% % if worm_frame_matrix(num_iteration,4)==0
% %     continue
% % end

worm_of_interest=worm_frame_matrix(num_iteration,2);
sprintf('Plotting video for worm %d out of %d worms', num_iteration, size(worm_frame_matrix,1))

str1 = "worm_";
str2 = num2str(worm_of_interest);
str3 = ".mat";

worm_image_data = append(str1,str2,str3);
load(worm_image_data)

%%%%%%%%%%%% plotting individual worm videos %%%%%%%%%%%%%%
all_corrected_projector_frames=loaded_variable.processed_decoded_camera_frames(current_tracks(worm_of_interest).Frames); %%%% finding the correct corresponding projector frame for the worm of interest

initial_frame=worm_frame_matrix(num_iteration,3);

if isempty(stimulus_intensities)
    continue
end

video_folder_path=fullfile(main_folder,video_folder_name);
    if ~exist(video_folder_path, 'dir')
        mkdir(fullfile(video_folder_path))
    end
    
    only_video_file_name="worm_"+num2str(worm_of_interest)+"_frame_"+num2str(initial_frame)+".avi";
    video_file_name = fullfile(video_folder_path,only_video_file_name);
    outputVideo = VideoWriter(video_file_name,'Uncompressed AVI');
    outputVideo.FrameRate = output_video_frame_rate;
    open(outputVideo);

        
for in_track_index=[(initial_frame-500):5:(initial_frame-100) (initial_frame-100):5:(initial_frame+100) (initial_frame+100):5:(initial_frame+500)]
    
    if (initial_frame-500)<=0 %%%% skip the iteration if the stim happens at the beginning of the video
        continue
    end
    
% % %     if (initial_frame+current_rails_dur)>size(current_tracks(worm_of_interest).Frames,1) %%%% skip the iteration if the stim happens at the end of the video
% % %         continue
% % %     end
    
    %%% being extra cautious; we are not plotting any videos when the worm is not tracked for 500 frames after stim onset
    if (initial_frame+500)>size(current_tracks(worm_of_interest).Frames,1) %%%% skip the iteration if the stim happens at the end of the video
        continue
    end
    
    if in_track_index>size(current_tracks(worm_of_interest).Frames,1)
        continue
    end
    
    current_correct_projector_frame=all_corrected_projector_frames(1,in_track_index);
    str11 = "Frame_";
    str12 = sprintf('%06d', current_correct_projector_frame);
    str13 = ".png";
    appended_str_projector = append(str11,str12,str13);
         
    curProjImage=imread(fullfile(main_folder,'ConvertedProjectorFrames',appended_str_projector'));

    centroid_x = double(round(current_tracks(worm_of_interest).Path(in_track_index,1)));
    centroid_y = double(round(current_tracks(worm_of_interest).Path(in_track_index,2)));
    image_top_left_corner_x = centroid_x-image_size(1)/2;
    image_top_left_corner_y = centroid_y-image_size(2)/2;
    image_bottom_right_corner_x = image_top_left_corner_x+image_size(1);
    image_bottom_right_corner_y = image_top_left_corner_y+image_size(2);
    cropped_image = imcrop(curProjImage, [image_top_left_corner_x, image_top_left_corner_y, (image_size-1)]);
    
    current_ellipse_ratio=current_tracks(worm_of_interest).EllipseRatio(in_track_index,1);
    all_ellipse_ratio=current_tracks(worm_of_interest).EllipseRatio(initial_frame-500:initial_frame+500,1);
    all_velocity=current_tracks(worm_of_interest).Velocity(initial_frame-500:initial_frame+500,1);
    current_velocity=current_tracks(worm_of_interest).Velocity(in_track_index,1);
    all_behavior=current_tracks(worm_of_interest).BehavioralAnnotation(1,:);
    
    path_of_worm=current_tracks(worm_of_interest).Path;
    stimulus_delivered=current_tracks(worm_of_interest).AlignedStimulus(in_track_index,10);
    
    %%%pad the image if necessary
    if image_top_left_corner_x < 1 || image_top_left_corner_y < 1
        %%%pad the front
        cropped_image = padarray(cropped_image, [max(1-image_top_left_corner_y,0), max(1-image_top_left_corner_x,0), 0], 0, 'pre');
    end
    if image_bottom_right_corner_x > size(curProjImage,2) || image_bottom_right_corner_y > size(curProjImage,1)
        %%%pad the end
        cropped_image = padarray(cropped_image, [max(image_bottom_right_corner_y-size(curProjImage,1)-1,0), max(image_bottom_right_corner_x-size(curProjImage,2)-1,0), 0], 0, 'post');
    end
    I = squeeze(worm_images(:,:,in_track_index));
    current_behavior = current_tracks(worm_of_interest).BehavioralAnnotation(in_track_index);
    if current_behavior < 1 || current_behavior > length(behavior_colors)
        current_behavior_color = [1 1 1];
    else
        current_behavior_color = behavior_colors(current_tracks(worm_of_interest).BehavioralAnnotation(in_track_index),:);
    end
    %%%%% convert the stimulus in uW/mm^2 to a color
    current_stimulus_along_centerline = current_tracks(worm_of_interest).AlignedStimulus(in_track_index,:);
    for centerline_point_index = 1:length(current_stimulus_along_centerline) %convert power intensity back to pixel intensity
        if current_stimulus_along_centerline(centerline_point_index) < 0
            %%% blue
            current_stimulus_along_centerline(centerline_point_index) = max(current_stimulus_along_centerline(centerline_point_index) / parameters.avgPowerBlue * 255, -255);
        else
            current_stimulus_along_centerline(centerline_point_index) = min(current_stimulus_along_centerline(centerline_point_index) / parameters.avgPowerRed * 255, 255);
        end
    end
    
% % %     %%% orignial code
% % %     plot_worm_frame_post_analysis_showing_path(I, squeeze(current_tracks(worm_of_interest).Centerlines(:,:,in_track_index)), current_behavior_color, ...
% % %             current_stimulus_along_centerline, cropped_image,in_track_index,current_ellipse_ratio);
    
    %%% code to display path in individual worm videos
    smoothed_path_x=smooth(current_tracks(worm_of_interest).Path(:,1),10);  %%% smoothing path with window 10
    smoothed_path_y=smooth(current_tracks(worm_of_interest).Path(:,2),10);
    
% % % % %     %%% to avoid any error if size of I do not match the size of debugimage
% % % % %     if size(I,1)~=size(cropped_image,1)
% % % % %         continue
% % % % %     end
    
    %%%% this works perfectly fine
% %     plot_worm_frame_post_analysis_showing_path(I, squeeze(current_tracks(worm_of_interest).Centerlines(:,:,in_track_index)), current_behavior_color, ...
% %             current_stimulus_along_centerline, cropped_image,in_track_index,current_ellipse_ratio,smoothed_path_x,smoothed_path_y,stimulus_delivered,path_history,demo_mode);
    
    %%% TRYING SUBPLOT
    plot_worm_frame_post_analysis_showing_path_subplot(I, squeeze(current_tracks(worm_of_interest).Centerlines(:,:,in_track_index)), current_behavior_color, ...
            stim_color, cropped_image,in_track_index,current_ellipse_ratio,smoothed_path_x,smoothed_path_y,stimulus_delivered,path_history,demo_mode,initial_frame,all_ellipse_ratio,worm_of_interest,current_rails_dur,all_velocity,current_velocity);
    
    hold on;
    hold off;
    writeVideo(outputVideo, getframe(gcf));
end

    close(outputVideo);
    video_file_name=convertStringsToChars(video_file_name);
    video_transcode(video_file_name);
    clc
end

end