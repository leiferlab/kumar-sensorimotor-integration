%%% To run this code you require workspace data. You can generate by
%%% running custom_optotap_stimulus_script or any of its variant. You need to make
%%% sure that when you run that file you mention these variables in 
%%% relevant_track_fields = {'BehavioralTransition','Frames',
%%% 'AlignedStimulus','EllipseRatio',%%% 'Path','Centerlines','UncertainTips'};
%%% It will generate video for "worm_of_interest". You can get the worm ID
%%% from processed.mkv file. Use 'analyze_particular_behavior'to study the
%%% behavior before and after the light stim. It will now analyze whole
%%% day data. 
%%% Last updated on 03/23/2021

clear
clc
close all

main_folder=('/projects/LEIFER/Sandeep/APIData/20220222_RunStimulateReversingWorms_Sandeep_AML496_lim4ACR2_10ulRet_blue/Data20220222_121519'); %% main folder path
load(fullfile(main_folder,'workspace_data_20220222_121519_AML496_180uW_10sec.mat')); %% save the workspace data inside the main folder

main_folder=('/projects/LEIFER/Sandeep/APIData/20220222_RunStimulateReversingWorms_Sandeep_AML496_lim4ACR2_10ulRet_blue/Data20220222_121519'); %% main folder path

%%%%%%%%%%%%%%%%%%%%%%%%%%%% user inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save_video=1;
analyze_indivial_worm=11;
analyze_whole_folder=1;
analyze_whole_day_data=11;
output_video_frame_rate=10;
behavior_to_study=2;
analyze_particular_behavior=11;
save_excel_file=11;
correct_folder_ID=11;
path_history=10*parameters.SampleRate;   %%%% how long the previous path should be drawn (in frames)
demo_mode=11;
%%

% % % % if analyze_indivial_worm==1

% % % % addpath(fullfile(main_folder,'/individual_worm_imgs'));
% % % % addpath(fullfile(main_folder,'/ConvertedProjectorFrames/'));
% % % % 
% % % % Centerlines=load(fullfile(main_folder,'/analysis/Centerlines.mat'));
% % % % Path=load(fullfile(main_folder,'/analysis/Path.mat'));
% % % % 
% % % % loaded_variable=load(fullfile(main_folder,'timestamps.mat'));
% % % % processed_decoded_camera_frames = loaded_variable.processed_decoded_camera_frames;

% % % % 
% % % % worm_of_interest=3; %%% e.g. 43 exampe: 141 for 20200702_193349 or 391 for 20200629_223802
% % % % 
% % % % str1 = "worm_";
% % % % str2 = num2str(worm_of_interest);
% % % % str3 = ".mat";
% % % % 
% % % % worm_image_data = append(str1,str2,str3);
% % % % load(worm_image_data)
% % % % 
% % % % %%%%%%%%% detecting stimulation when worm is turning %%%%%%%%%%%%%%%%%%%%%%
% % % % data_stim=find(current_tracks(worm_of_interest).AlignedStimulus(:,10)/max(current_tracks(worm_of_interest).AlignedStimulus(:,10))*2==2); %% when LED were ON
% % % % data_turn=find(current_tracks(worm_of_interest).BehavioralAnnotation==behavior_to_study); %% when worm was turning
% % % % correct_stims=intersect(data_stim,data_turn);
% % % % 
% % % % correct_stim_final=zeros(1,size(worm_images,3));
% % % % correct_stim_final(correct_stims)=2;
% % % % 
% % % % data_turn=find(current_tracks(worm_of_interest).BehavioralAnnotation==behavior_to_study); %% when worm is turning
% % % % data_turn_final=zeros(1,size(current_tracks(worm_of_interest).AlignedStimulus,1));
% % % % data_turn_final(data_turn)=2; %%% this is the correct stimulation array
% % % % 
% % % % data_reverse=find(current_tracks(worm_of_interest).BehavioralAnnotation==3); %% when worm is reversing
% % % % data_reverse_final=zeros(1,size(current_tracks(worm_of_interest).AlignedStimulus,1));
% % % % data_reverse_final(data_reverse)=3; %%% this is the correct stimulation array
% % % % 
% % % % figure; 
% % % % plot(current_tracks(worm_of_interest).AlignedStimulus(:,1)/current_stim_power*4,'-r','DisplayName','Opto. stim.'); 
% % % % hold on; 
% % % % % plot(current_tracks(worm_of_interest).EllipseRatio,'-b','DisplayName','Ellipse ratio')
% % % % % hold on;
% % % % % plot(data_turn_final,'-g','DisplayName','Turning state')
% % % % % hold on;
% % % % % plot(data_reverse_final,'-k','DisplayName','Reversal state')
% % % % % hold on;
% % % % % plot(current_tracks(worm_of_interest).Velocity*10,'-c','DisplayName','Velocity'); hold on;
% % % % % plot(current_tracks(worm_of_interest).Direction/20,'-k','DisplayName','Direction'); hold on;
% % % % % plot(current_tracks(worm_of_interest).BehavioralAnnotation,'-g','DisplayName','Beh. state'); hold on;
% % % % % plot(correct_stim_final,'--k','DisplayName','correct stim.'); % hold on;
% % % % % yline(parameters.EllipseRatioTriggerPoint,'--m','Linewidth',2,'DisplayName','E.R. trigger point');
% % % % ylim([0 10])
% % % % % xlim([0 3000])
% % % % xlabel('Frame count') % x-axis label
% % % % ylabel('Ellipse ratio') % y-axis label
% % % % title(['Worm index = ', num2str(worm_of_interest)])
% % % % ax = gca;
% % % % ax.FontSize = 13;
% % % % legend('show','Location','southeast','FontSize',12);
% % % % box off; axis tight
% % % % 
% % % % %%
% % % % pause(1)
% % % % %%%%%%%%%%%% plotting individual worm videos %%%%%%%%%%%%%%
% % % % all_corrected_projector_frames=loaded_variable.processed_decoded_camera_frames(current_tracks(worm_of_interest).Frames); %%%% finding the correct corresponding projector frame for the worm of interest
% % % % 
% % % % ImageSize = 2*round(parameters.TrackingWormBoxHalfSize * parameters.CameraPixeltommConversion);
% % % % image_size = [ImageSize, ImageSize];
% % % % 
% % % % initial_frame=509-500;
% % % % 
% % % % figure;
% % % % if save_video==1
% % % %     if ~exist([main_folder, filesep, 'individual_worm_videos/turning_180uW_60uW_10sec'], 'dir')
% % % %         mkdir([main_folder, filesep, 'individual_worm_videos/turning_180uW_60uW_10sec'])
% % % %     end
% % % %     video_file_name = fullfile([main_folder, filesep, 'individual_worm_videos/turning_180uW_60uW_10sec', filesep, 'worm_', num2str(worm_of_interest), '_frame_',num2str(initial_frame),'.avi']);
% % % %     outputVideo = VideoWriter(video_file_name,'Uncompressed AVI');
% % % %     outputVideo.FrameRate = output_video_frame_rate;
% % % %     open(outputVideo);
% % % % end
% % % %         
% % % % for in_track_index=initial_frame:5:initial_frame+1000
% % % %     in_track_index
% % % %     current_correct_projector_frame=all_corrected_projector_frames(in_track_index);
% % % %     str11 = "Frame_";
% % % %     str12 = sprintf('%06d', current_correct_projector_frame);
% % % %     str13 = ".png";
% % % %     appended_str_projector = append(str11,str12,str13);
% % % %          
% % % %     curProjImage=imread(fullfile(main_folder,'ConvertedProjectorFrames',appended_str_projector'));
% % % % 
% % % %     centroid_x = double(round(current_tracks(worm_of_interest).Path(in_track_index,1)));
% % % %     centroid_y = double(round(current_tracks(worm_of_interest).Path(in_track_index,2)));
% % % %     image_top_left_corner_x = centroid_x-image_size(1)/2;
% % % %     image_top_left_corner_y = centroid_y-image_size(2)/2;
% % % %     image_bottom_right_corner_x = image_top_left_corner_x+image_size(1);
% % % %     image_bottom_right_corner_y = image_top_left_corner_y+image_size(2);
% % % %     cropped_image = imcrop(curProjImage, [image_top_left_corner_x, image_top_left_corner_y, (image_size-1)]);
% % % %     
% % % %     current_ellipse_ratio=current_tracks(worm_of_interest).EllipseRatio(in_track_index,1);
% % % %     %pad the image if necessary
% % % %     if image_top_left_corner_x < 1 || image_top_left_corner_y < 1
% % % %         %pad the front
% % % %         cropped_image = padarray(cropped_image, [max(1-image_top_left_corner_y,0), max(1-image_top_left_corner_x,0), 0], 0, 'pre');
% % % %     end
% % % %     if image_bottom_right_corner_x > size(curProjImage,2) || image_bottom_right_corner_y > size(curProjImage,1)
% % % %         %pad the end
% % % %         cropped_image = padarray(cropped_image, [max(image_bottom_right_corner_y-size(curProjImage,1)-1,0), max(image_bottom_right_corner_x-size(curProjImage,2)-1,0), 0], 0, 'post');
% % % %     end
% % % %     I = squeeze(worm_images(:,:,in_track_index));
% % % %     current_behavior = current_tracks(worm_of_interest).BehavioralAnnotation(in_track_index);
% % % %     if current_behavior < 1 || current_behavior > length(behavior_colors)
% % % %         current_behavior_color = [1 1 1];
% % % %     else
% % % %         current_behavior_color = behavior_colors(current_tracks(worm_of_interest).BehavioralAnnotation(in_track_index),:);
% % % %     end
% % % %     %%%%% convert the stimulus in uW/mm^2 to a color
% % % %     current_stimulus_along_centerline = current_tracks(worm_of_interest).AlignedStimulus(in_track_index,:);
% % % %     for centerline_point_index = 1:length(current_stimulus_along_centerline) %convert power intensity back to pixel intensity
% % % %         if current_stimulus_along_centerline(centerline_point_index) < 0
% % % %             % blue
% % % %             current_stimulus_along_centerline(centerline_point_index) = max(current_stimulus_along_centerline(centerline_point_index) / parameters.avgPowerBlue * 255, -255);
% % % %         else
% % % %             current_stimulus_along_centerline(centerline_point_index) = min(current_stimulus_along_centerline(centerline_point_index) / parameters.avgPowerRed * 255, 255);
% % % %         end
% % % %     end
% % % %     
% % % %     plot_worm_frame_post_analysis(I, squeeze(current_tracks(worm_of_interest).Centerlines(:,:,in_track_index)), current_behavior_color, ... 
% % % %             current_stimulus_along_centerline, cropped_image,in_track_index,current_ellipse_ratio);
% % % %     
% % % % %     plot_worm_frame_post_analysis(I, squeeze(current_tracks(worm_of_interest).Centerlines(:,:,in_track_index)), current_behavior_color, ...
% % % % %         current_tracks(worm_of_interest).UncertainTips(in_track_index), current_stimulus_along_centerline, ...
% % % % %         cropped_image,in_track_index,current_ellipse_ratio);
% % % % %     
% % % %     hold on;
% % % %     hold off;
% % % %     pause(0.01);
% % % %     if save_video==1
% % % %     writeVideo(outputVideo, getframe(gcf));
% % % %     end
% % % % end
% % % % 
% % % % if save_video==1
% % % %     close(outputVideo);
% % % %     video_transcode(video_file_name);
% % % % end
% % % % end
% % % % 
%%
% h=figure
close all
figure1=figure('Position', [40 40 1120 640]);

if analyze_whole_folder==1

addpath(fullfile(main_folder,'/individual_worm_imgs'));
addpath(fullfile(main_folder,'/ConvertedProjectorFrames/'));

% % Centerlines=load(fullfile(main_folder,'/analysis/Centerlines.mat'));
% % Path=load(fullfile(main_folder,'/analysis/Path.mat'));

loaded_variable=load(fullfile(main_folder,'timestamps.mat'));
processed_decoded_camera_frames = loaded_variable.processed_decoded_camera_frames;

if save_excel_file==1
% row_id=[1:1:size(stim_while_turning_peaks_array,1)]';
% stim_while_turning_peaks_array=[row_id stim_while_turning_peaks_array];
stim_while_turning_peaks_array=[stim_while_turning_peaks_array];
indices_with_zero_rows=find(~all(stim_while_turning_peaks_array(:,2:end)==0,2));
stim_while_turning_peaks_array_excel_data=stim_while_turning_peaks_array(indices_with_zero_rows,:);
excel_filename = 'new_AKS449.3.1.dxAML67_data_20210207_160238_180uW_10sec.xlsx';
writematrix(stim_while_turning_peaks_array_excel_data,excel_filename)
end

ImageSize = 2*round(parameters.TrackingWormBoxHalfSize * parameters.CameraPixeltommConversion);
image_size = [ImageSize, ImageSize];

if contains(main_folder, 'Full')
    testtype='rails';
    worm_frame_matrix=count_number_of_reversal_trial_rails;
elseif contains(main_folder,'Turn')
    testtype='turns';
    worm_frame_matrix=count_number_of_reversal_trial_turns;
elseif contains(main_folder,'Reversing')
    testtype='turns';
    worm_frame_matrix=count_number_of_reversal_trial_turns;
end

%%% in case you want to plot data from another variable
% load('/projects/LEIFER/Sandeep/APIData/20210702_RunRailsTriggeredByTurning_Sandeep_AKS_483.7.e_mec4_Chrimson_10ulRet_red/Data20210702_170102/plot_additional_individual_worm_videos.mat')
% worm_frame_matrix=plot_additional_individual_worm_videos;
% worm_frame_matrix=worm_frame_matrix(47,:);
% worm_frame_matrix=[worm_frame_matrix 
%     1 195 1994 1];

% worm_frame_matrix=worm_frame_matrix([15 36 37],:);
% worm_frame_matrix=[1 1 1 7502; 1 3 1 3842; 1 16 1 7760];
% % for i=1:size(current_tracks,2)
% %     worm_frame_matrix(i,1)=1;
% %     worm_frame_matrix(i,2)=i;
% %     worm_frame_matrix(i,3)=1;
% %     worm_frame_matrix(i,4)=size(current_tracks(i).Frames,1);
% % end

worm_frame_matrix=worm_frame_matrix([4 5 12 13 17 20 29 31 38 70],:); %% 4 5 12 13 17 20 29 31 38 70

%%
for num_iteration=1:size(worm_frame_matrix,1) 
% % for num_iteration=106 %%%%% for debugging

% % %%%%%% save video only when certain criteria is met
% % if worm_frame_matrix(num_iteration,4)==0
% %     continue
% % end

worm_of_interest=worm_frame_matrix(num_iteration,2);
sprintf('Analyzing worm %d out of %d worms', num_iteration, size(worm_frame_matrix,1))

str1 = "worm_";
str2 = num2str(worm_of_interest);
str3 = ".mat";

worm_image_data = append(str1,str2,str3);
load(worm_image_data)

%%%%%%%%%%%% plotting individual worm videos %%%%%%%%%%%%%%
all_corrected_projector_frames=loaded_variable.processed_decoded_camera_frames(current_tracks(worm_of_interest).Frames); %%%% finding the correct corresponding projector frame for the worm of interest

initial_frame=worm_frame_matrix(num_iteration,3);
% % final_frame=worm_frame_matrix(num_iteration,4);

if save_video==1
    if ~exist([main_folder, filesep, 'individual_worm_videos/no_projector_image_version/rails_180UW_3s'], 'dir')
        mkdir([main_folder, filesep, 'individual_worm_videos/no_projector_image_version/rails_180UW_3s'])
    end
    video_file_name = fullfile([main_folder, filesep, 'individual_worm_videos/no_projector_image_version/rails_180UW_3s', filesep, 'worm_', num2str(worm_of_interest), '_frame_',num2str(initial_frame),'.avi']);
    outputVideo = VideoWriter(video_file_name,'Uncompressed AVI');
    outputVideo.FrameRate = output_video_frame_rate;
    open(outputVideo);
end
        
for in_track_index=(initial_frame-500):5:(initial_frame+500)
% for in_track_index=(initial_frame):5:(final_frame)
% for in_track_index=1:5:5081
    
    if (initial_frame-500)<=0 %%%% skip the iteration if the stim happens at the beginning of the video
        continue
    end
    
    if (initial_frame+current_rails_dur)>size(current_tracks(worm_of_interest).Frames,1) %%%% skip the iteration if the stim happens at the end of the video
        continue
    end
    
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
    
% % % %     %%%% this works perfectly fine
% % % %     plot_worm_frame_post_analysis_showing_path(I, squeeze(current_tracks(worm_of_interest).Centerlines(:,:,in_track_index)), current_behavior_color, ...
% % % %             current_stimulus_along_centerline, cropped_image,in_track_index,current_ellipse_ratio,smoothed_path_x,smoothed_path_y,stimulus_delivered,path_history,demo_mode);
    
    %%%% trying to not show projector image
%     plot_worm_frame_post_analysis_showing_path_no_projector_image(I, squeeze(current_tracks(worm_of_interest).Centerlines(:,:,in_track_index)), current_behavior_color, ...
%             current_stimulus_along_centerline,in_track_index,current_ellipse_ratio,smoothed_path_x,smoothed_path_y,stimulus_delivered,path_history,demo_mode);
   
    %%% TRYING SUBPLOT
% % %     current_rails_dur=90;
% % %     plot_worm_frame_post_analysis_showing_path_subplot(I, squeeze(current_tracks(worm_of_interest).Centerlines(:,:,in_track_index)), current_behavior_color, ...
% % %             stim_color, cropped_image,in_track_index,current_ellipse_ratio,smoothed_path_x,smoothed_path_y,stimulus_delivered,path_history,demo_mode,initial_frame,all_ellipse_ratio,worm_of_interest,current_rails_dur,all_velocity,current_velocity);
    
% % % %     plot_worm_frame_post_analysis(I, squeeze(current_tracks(worm_of_interest).Centerlines(:,:,in_track_index)), current_behavior_color, ...
% % % %         current_tracks(worm_of_interest).UncertainTips(in_track_index), current_stimulus_along_centerline, ...
% % % %         cropped_image,in_track_index,current_ellipse_ratio);
    
    hold on;
    hold off;
    if save_video==1
    writeVideo(outputVideo, getframe(gcf));
    end
end

if save_video==1
    close(outputVideo);
    video_transcode(video_file_name);
end
end
end

%%
%%
% % % % % % % if analyze_whole_day_data==1
% % % % % % % 
% % % % % % % load('reversal_duration_folder_20210207_449_31d_AML67_red_60uW_turn.mat')
% % % % % % % 
% % % % % % % index_for_correct_folder = any(reversal_duration_folder_20210207_449_31d_AML67_red_60uW_turn(:,1) == correct_folder_ID,2);
% % % % % % % reversal_duration_folder=reversal_duration_folder_20210207_449_31d_AML67_red_60uW_turn(index_for_correct_folder,:);
% % % % % % % 
% % % % % % % if save_excel_file==1
% % % % % % % excel_filename = 'new_AKS449.3.1.dxAML67_data_20210207_160238_60uW_10sec.xlsx';
% % % % % % % writematrix(stim_while_turning_peaks_array_excel_data,excel_filename)
% % % % % % % end
% % % % % % % 
% % % % % % % ImageSize = 2*round(parameters.TrackingWormBoxHalfSize * parameters.CameraPixeltommConversion);
% % % % % % % image_size = [ImageSize, ImageSize];
% % % % % % % 
% % % % % % % addpath(fullfile(main_folder,'/individual_worm_imgs'));
% % % % % % % addpath(fullfile(main_folder,'/ConvertedProjectorFrames/'));
% % % % % % % 
% % % % % % % Centerlines=load(fullfile(main_folder,'/analysis/Centerlines.mat'));
% % % % % % % Path=load(fullfile(main_folder,'/analysis/Path.mat'));
% % % % % % % 
% % % % % % % loaded_variable=load(fullfile(main_folder,'timestamps.mat'));
% % % % % % % processed_decoded_camera_frames = loaded_variable.processed_decoded_camera_frames;
% % % % % % % 
% % % % % % % for xx=1:size(reversal_duration_folder,1) 
% % % % % % %     
% % % % % % % sprintf('Analyzing worm %d out of %d worms', xx, size(reversal_duration_folder,1))
% % % % % % % 
% % % % % % % worm_of_interest=reversal_duration_folder(xx,2);
% % % % % % % 
% % % % % % % str1 = "worm_";
% % % % % % % str2 = num2str(worm_of_interest);
% % % % % % % str3 = ".mat";
% % % % % % % 
% % % % % % % worm_image_data = append(str1,str2,str3);
% % % % % % % load(worm_image_data)
% % % % % % % 
% % % % % % % %%%%%%%%%%%% plotting individual worm videos %%%%%%%%%%%%%%
% % % % % % % all_corrected_projector_frames=loaded_variable.processed_decoded_camera_frames(current_tracks(worm_of_interest).Frames); %%%% finding the correct corresponding projector frame for the worm of interest
% % % % % % % 
% % % % % % % peaks_data=reversal_duration_folder(xx,3);
% % % % % % % nonzero_peaks_data=nonzeros(peaks_data);
% % % % % % % 
% % % % % % % if size(nonzero_peaks_data,1)==0  %%% skip the iteration if there are no non-zero elemets
% % % % % % %     continue
% % % % % % % end
% % % % % % % 
% % % % % % % for mn=1:size(nonzero_peaks_data,1)
% % % % % % % initial_frame=nonzero_peaks_data(mn);
% % % % % % % 
% % % % % % % if save_video==1
% % % % % % %     if ~exist([main_folder, filesep, 'individual_worm_videos/turning_180uW_60uW_10sec'], 'dir')
% % % % % % %         mkdir([main_folder, filesep, 'individual_worm_videos/turning_180uW_60uW_10sec'])
% % % % % % %     end
% % % % % % %     video_file_name = fullfile([main_folder, filesep, 'individual_worm_videos/turning_180uW_60uW_10sec', filesep, 'worm_', num2str(worm_of_interest), '_frame_',num2str(initial_frame),'.avi']);
% % % % % % %     outputVideo = VideoWriter(video_file_name,'Uncompressed AVI');
% % % % % % %     outputVideo.FrameRate = output_video_frame_rate;
% % % % % % %     open(outputVideo);
% % % % % % % end
% % % % % % %         
% % % % % % % for in_track_index=(initial_frame-500):5:(initial_frame+500)
% % % % % % %     
% % % % % % %     if (initial_frame-500)<=0 %%%% skip the iteration if the stim happens at the beginning of the video
% % % % % % %         disp('Initial_frame<500')
% % % % % % %         continue
% % % % % % %     end
% % % % % % %     
% % % % % % %     if (initial_frame+500)>size(current_tracks(worm_of_interest).Frames,1) %%%% skip the iteration if the stim happens at the end of the video
% % % % % % %         disp('Final_frame>500')
% % % % % % %         continue
% % % % % % %     end
% % % % % % %     
% % % % % % %     current_correct_projector_frame=all_corrected_projector_frames(in_track_index);
% % % % % % %     str11 = "Frame_";
% % % % % % %     str12 = sprintf('%06d', current_correct_projector_frame);
% % % % % % %     str13 = ".png";
% % % % % % %     appended_str_projector = append(str11,str12,str13);
% % % % % % %          
% % % % % % %     curProjImage=imread(fullfile(main_folder,'ConvertedProjectorFrames',appended_str_projector'));
% % % % % % % 
% % % % % % %     centroid_x = double(round(current_tracks(worm_of_interest).Path(in_track_index,1)));
% % % % % % %     centroid_y = double(round(current_tracks(worm_of_interest).Path(in_track_index,2)));
% % % % % % %     image_top_left_corner_x = centroid_x-image_size(1)/2;
% % % % % % %     image_top_left_corner_y = centroid_y-image_size(2)/2;
% % % % % % %     image_bottom_right_corner_x = image_top_left_corner_x+image_size(1);
% % % % % % %     image_bottom_right_corner_y = image_top_left_corner_y+image_size(2);
% % % % % % %     cropped_image = imcrop(curProjImage, [image_top_left_corner_x, image_top_left_corner_y, (image_size-1)]);
% % % % % % %     
% % % % % % %     current_ellipse_ratio=current_tracks(worm_of_interest).EllipseRatio(in_track_index,1);
% % % % % % %     %%%pad the image if necessary
% % % % % % %     if image_top_left_corner_x < 1 || image_top_left_corner_y < 1
% % % % % % %         %%%pad the front
% % % % % % %         cropped_image = padarray(cropped_image, [max(1-image_top_left_corner_y,0), max(1-image_top_left_corner_x,0), 0], 0, 'pre');
% % % % % % %     end
% % % % % % %     if image_bottom_right_corner_x > size(curProjImage,2) || image_bottom_right_corner_y > size(curProjImage,1)
% % % % % % %         %%%pad the end
% % % % % % %         cropped_image = padarray(cropped_image, [max(image_bottom_right_corner_y-size(curProjImage,1)-1,0), max(image_bottom_right_corner_x-size(curProjImage,2)-1,0), 0], 0, 'post');
% % % % % % %     end
% % % % % % %     I = squeeze(worm_images(:,:,in_track_index));
% % % % % % %     current_behavior = current_tracks(worm_of_interest).BehavioralAnnotation(in_track_index);
% % % % % % %     if current_behavior < 1 || current_behavior > length(behavior_colors)
% % % % % % %         current_behavior_color = [1 1 1];
% % % % % % %     else
% % % % % % %         current_behavior_color = behavior_colors(current_tracks(worm_of_interest).BehavioralAnnotation(in_track_index),:);
% % % % % % %     end
% % % % % % %     %%%%% convert the stimulus in uW/mm^2 to a color
% % % % % % %     current_stimulus_along_centerline = current_tracks(worm_of_interest).AlignedStimulus(in_track_index,:);
% % % % % % %     for centerline_point_index = 1:length(current_stimulus_along_centerline) %convert power intensity back to pixel intensity
% % % % % % %         if current_stimulus_along_centerline(centerline_point_index) < 0
% % % % % % %             %%% blue
% % % % % % %             current_stimulus_along_centerline(centerline_point_index) = max(current_stimulus_along_centerline(centerline_point_index) / parameters.avgPowerBlue * 255, -255);
% % % % % % %         else
% % % % % % %             current_stimulus_along_centerline(centerline_point_index) = min(current_stimulus_along_centerline(centerline_point_index) / parameters.avgPowerRed * 255, 255);
% % % % % % %         end
% % % % % % %     end
% % % % % % %     
% % % % % % %         plot_worm_frame_post_analysis(I, squeeze(current_tracks(worm_of_interest).Centerlines(:,:,in_track_index)), current_behavior_color, ... 
% % % % % % %             current_stimulus_along_centerline, cropped_image,in_track_index,current_ellipse_ratio);
% % % % % % %     
% % % % % % % % % % %     plot_worm_frame_post_analysis(I, squeeze(current_tracks(worm_of_interest).Centerlines(:,:,in_track_index)), current_behavior_color, ...
% % % % % % % % % % %         current_tracks(worm_of_interest).UncertainTips(in_track_index), current_stimulus_along_centerline, ...
% % % % % % % % % % %         cropped_image,in_track_index,current_ellipse_ratio);
% % % % % % %     
% % % % % % %     hold on;
% % % % % % %     hold off;
% % % % % % %     if save_video==1
% % % % % % %     writeVideo(outputVideo, getframe(gcf));
% % % % % % %     end
% % % % % % % end
% % % % % % % 
% % % % % % % if save_video==1
% % % % % % %     close(outputVideo);
% % % % % % %     video_transcode(video_file_name);
% % % % % % % end
% % % % % % % % % % clearvars I squeeze(current_tracks(worm_of_interest).Centerlines(:,:,in_track_index)) current_behavior_color ... 
% % % % % % % % % %             current_stimulus_along_centerline cropped_image in_track_index current_ellipse_ratio;
% % % % % % % end
% % % % % % % end
% % % % % % % % % % end
% % % % % % % end
% % % % % % % %%
% % % % % % % if analyze_particular_behavior==1
% % % % % % % 
% % % % % % % for worm_of_interest=1:size(stim_while_turning_peaks_array,1) 
% % % % % % % % for worm_of_interest=1:5  %%%%% for debugging
% % % % % % %     
% % % % % % % peaks_data=stim_while_turning_peaks_array(worm_of_interest,:);
% % % % % % % nonzero_peaks_data=nonzeros(peaks_data);
% % % % % % % 
% % % % % % % if size(nonzero_peaks_data,1)==0  %%% skip the iteration if there are no non-zero elemets
% % % % % % %     continue
% % % % % % % end
% % % % % % % 
% % % % % % % for mn=1:size(nonzero_peaks_data,1)
% % % % % % % % for mn=1
% % % % % % % initial_frame=nonzero_peaks_data(mn);
% % % % % % % 
% % % % % % %  current_velocity_beh_all_tracks=[];
% % % % % % %  counter=0;
% % % % % % % for in_track_index=(initial_frame-500):1:(initial_frame+500)
% % % % % % %     counter=counter+1;
% % % % % % %     if (initial_frame-500)<0 %%%% skip the iteration if the stim happens at the beginning of the video
% % % % % % %         continue
% % % % % % %     end
% % % % % % %     
% % % % % % %     if (initial_frame+500)>size(current_tracks(worm_of_interest).Frames,1) %%%% skip the iteration if the stim happens at the end of the video
% % % % % % %         continue
% % % % % % %     end
% % % % % % %     
% % % % % % %     current_velocity_beh_track(counter,mn)=current_tracks(worm_of_interest).Velocity(in_track_index,1);
% % % % % % %     
% % % % % % % end
% % % % % % % 
% % % % % % % mean_velocity_beh_worm(:,worm_of_interest)=mean(current_velocity_beh_track,2);
% % % % % % % 
% % % % % % % end
% % % % % % % end
% % % % % % % 
% % % % % % % mean_velocity_beh_worm( :, ~any(mean_velocity_beh_worm,1) ) = [];  %%%deleting zero columns
% % % % % % % mean_velocity_all_worms=mean(mean_velocity_beh_worm,2);
% % % % % % % std_velocity_all_worms=std(mean_velocity_beh_worm,0,2);
% % % % % % % x_axis=[1:1:size(mean_velocity_all_worms,1)];
% % % % % % % figure;
% % % % % % % errorbar(x_axis,mean_velocity_all_worms,std_velocity_all_worms)
% % % % % % % hold on;
% % % % % % % plot(x_axis,mean_velocity_all_worms,'-r','LineWidth',2)
% % % % % % % 
% % % % % % % mean_ellipse_ratio_all_worms_AML17=mean_velocity_all_worms;
% % % % % % % std_ellipse_ratio_all_worms_AML17=std_velocity_all_worms;
% % % % % % % end
% % % % % % % 
% % % % % % % % % % % %%
% % % % % % % % % % % %%%% plotting elipse ratio
% % % % % % % % % % % allratios = vertcat(current_tracks.EllipseRatio);
% % % % % % % % % % % figure;
% % % % % % % % % % % histogram(allratios,20,'facecolor','c');
% % % % % % % % % % % hold on;
% % % % % % % % % % % xline(parameters.EllipseRatioTriggerPoint,'--k','Linewidth',2);
% % % % % % % % % % % ax = gca;
% % % % % % % % % % % ax.FontSize = 13;
% % % % % % % % % % % xlabel('Ellipse ratio') % x-axis label
% % % % % % % % % % % ylabel('Number of occurences') % y-axis label
% % % % % % % % % % % title('Distribution of ellipse ratio');
% % % % % % % % % % % axis tight;
% % % % % % % % % % % box off;
% % % % % % % % % % % %%
% % % % % % % % % % % 
% % % % % % % % % % % %%%%%%%%% to plot all the peak widths
% % % % % % % % % % % 
% % % % % % % % % % % for worm_id=2:size(current_tracks,2)
% % % % % % % % % % %     mid_cline_index=size(current_tracks(worm_id).AlignedStimulus,2)/2;
% % % % % % % % % % %     LED_data=current_tracks(worm_id).AlignedStimulus(:,mid_cline_index);
% % % % % % % % % % %     [pks, locs,width,~]=findpeaks(LED_data,'Annotate','extents');
% % % % % % % % % % %     
% % % % % % % % % % %     if( isempty(width) )
% % % % % % % % % % %         width=NaN;
% % % % % % % % % % %     end
% % % % % % % % % % % 
% % % % % % % % % % %     width_cell{worm_id,:}=width;
% % % % % % % % % % % end
% % % % % % % % % % % 
% % % % % % % % % % % width_array=cell2mat(width_cell);
% % % % % % % % % % % figure;
% % % % % % % % % % % histogram(width_array,30);
% % % % % % % % % % % xlabel('LED stim width') % x-axis label
% % % % % % % % % % % ylabel('Number of occurences') % y-axis label
% % % % % % % % % % % title('Distribution of LED stim width');
% % % % % % % % % % % axis tight;
% % % % % % % % % % % box off;
% % % % % % % % % % % ax = gca;
% % % % % % % % % % % ax.FontSize = 13;