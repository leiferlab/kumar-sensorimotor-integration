%%%% This code is used in the analysis for generating data for table 1 and
%%%% figure 4 for AML67 and AML470. The number of instances in which reversals
%%%% are generated for closed loop and open loop assays are saved in th evariables
%%%% count_number_of_reversal_folder_turns and
%%%% count_number_of_reversal_folder_rails respectively.
%%%% Updated: 08/24/2021

addpath(genpath('/projects/LEIFER/Sandeep/github/ProjectAPI'));

close all
clear
clc

%%%% run for the last folder in case of examining just one folder for
%%%% dataset of 20220120
% % main_folder=('/projects/LEIFER/Sandeep/APIData/20220222_RunStimulateReversingWorms_Sandeep_AML496_lim4ACR2_10ulRet_blue'); %%% for AML496 data
main_folder=('/projects/LEIFER/Sandeep/APIData/20220120_RunStimulateReversingWorms_Sandeep_AML499_lim4ACR2_mec4Chrimson_10ulRet_blue'); %%% for AML499 data

cd(main_folder)

%%%%% master input from the user
random_seed_velocity_heatmap=1;
behavior_sort_threshold=3; %% 1=for, 2=turn, 3=rev
behavior_comparision_mode=2; %% 1:<=, 2:==, 3: >=
velocity_sort_threshold=0;
velocity_comparision_mode=1; %% 1:<=, 2:==, 3: >=
analysis_new_data=1;
analysis_workspace_data=11;

%%%%%%%%%%%%%%%%%%
if analysis_workspace_data==1
%    load('/projects/LEIFER/Sandeep/Publications/Kumar_etal_2022/workspace_datasets/workspace_data_AML67_20200902_stim_on_turns.mat') 
    load('/projects/LEIFER/Sandeep/Publications/Kumar_etal_2022/workspace_datasets/workspace_data_AML67_20200902_stim_on_fwd.mat') 
    analysis_new_data=11;
    analysis_workspace_data=1;
end

%%%%%%% user input parameters %%%%%%
test_stimulus_duration=10; %%% if you want to test a specific stim duration
stim_threshold_min=0.01; %%%% min stim power 
stim_threshold_max=10;   %%%% max stim power
plotting_graphs=1;       %% to avoid plotting graphs
save_mat_file=11;          %%%% to save the .mat file
plot_all_beh_ratios=1;
plot_various_int=11;
plot_comparison_reverse_beh=11;
plot_behavior_heatmaps=11;
plot_velocity_traces=11;
plot_ellipse_ratio_traces=11;
plot_fwd_to_turn_heatmaps=11;
plot_reverse_to_turn_heatmaps=11;
count_prob_reversal=11;
number_of_rastors=30;
plot_velocity_heatmaps=11;
plot_ellipse_ratio_heatmaps=11;
random_seed_open_loop_heatmap=14; %%% use 22 for open loop
random_seed_closed_loop_heatmap=39;
random_seed_reverse_to_turn_heatmap=3;
random_seed_forward_to_turn_heatmap=15;

excel_filename = 'AKS_483.7.e_20210723_175546_turns_20uW_3sec.csv';

%%%%%%% optional parameters to help with debugging
plot_head_tail_switch_instances=11;

%%%%%%% determine stim color %%%%%%%%
if contains(main_folder, 'blue')
    stim_color=[0 0 1];
elseif contains(main_folder,'purple')
    stim_color=[1 0 1];
elseif contains(main_folder,'red')
    stim_color=[1 0 0]; 
else
    stim_color=[0 0 1];
end

%%%%%%%%%%%%%%%%% Input parameters %%%%%%%%%%%%%%%%
behavior_to_study=3;  %%%  2= turn, 3=reverse 
two_stim_delivered=11; %%% when both blue and red stim are delivered simultaneously
opto_stim_color=3; %%% red=1, green=2, blue=3

max_allowed_stim_duration_diff=10; %%% this parameter is used for sorting stim of diff widths   
minimum_length_of_worm=0.55;  %%%% 
ellipse_ratio_threshold=3.1;  %%% use 3.6 as default value
negative_vel_threshold=-0.11;
frames_in_reversal_threshold=1; %%%% this is the time for which worms must reverse
test_low_ellipse_ratio_worm_threshold=3.9; %%% if the mean ellipse ratio of certian frames below this threshold then that trial is ignored
additional_seconds_for_visualization=2;   %%% In case you want to visualize data for longer than stim duration. Keep it default to 0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if contains(main_folder, 'Full')
    testtype='rails';
    analysis_rails=1;
    analysis_stim_while_turning=11;
elseif contains(main_folder,'Turn')
    testtype='turns';
    analysis_rails=11;
    analysis_stim_while_turning=1;
elseif contains(main_folder,'Reversing')
    testtype='turns';
    analysis_rails=11;
    analysis_stim_while_turning=1;
else
    analysis_rails=1;
    analysis_stim_while_turning=11;  
end

%%%% This if loop will pick the correct behavior to select whether it is reversing or turning.
if contains(main_folder, 'Reversing')
    correct_behavior_threshold=negative_vel_threshold;
    analysis_stim_while_turning=1;
    behavior_to_study=3;
elseif contains(main_folder,'Turn')
    correct_behavior_threshold=ellipse_ratio_threshold;
    analysis_stim_while_turning=1;
    behavior_to_study=2;
end

if contains(main_folder,'TwoColors')
    two_stim_delivered=1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if analysis_new_data==1
%%% ask for folders
load('reference_embedding.mat')

%%% if plotting individual videos then use this instead:
relevant_track_fields = {'BehavioralTransition','Frames','AlignedStimulus','EllipseRatio','Path','Centerlines','Velocity','Length','MeanAngle'};
% relevant_track_fields = {'BehavioralTransition','Frames','AlignedStimulus','Velocity','EllipseRatio','Length','Centerlines','MeanAngle'};

%%%% select folders
folders = getfoldersGUI();
% % % load('/projects/LEIFER/Sandeep/Publications/Liu_etal_2021/reviewer_comments/tracking_frame_rate/all_open_loop/folders.mat')

parameters = load_parameters(folders{end}); % load parameters for the first folder

%%
%%%%% constants
normalized_stimuli = 1; %delta function
fps = parameters.SampleRate;
stim_similarity_thresh = 1.4;  %%% use 1.1 % combine data in case the intensity is within some range (to take care of box to box variation)
number_of_behaviors = max(L(:)-1);
rails_durations=parameters.RailsDurations;
rails_durations=sort(rails_durations);
rails_durations=unique(rails_durations);  %%% in case, in th elabview program same time was entered twice

index_durations=find(rails_durations==test_stimulus_duration);
for rails_dur=index_durations   %%% to debug or analyze particular stim duration

% for rails_dur=1:size(rails_durations,2) %%% for analyzing all the stim durations

rails_dur_itr=rails_durations(1,rails_dur);
current_rails_dur=rails_dur_itr*parameters.SampleRate;

vel_xlim_turns_1=501; %%% defining the window in which I will look for negative vel for stim on turns
vel_xlim_turns_2=vel_xlim_turns_1+current_rails_dur; %%% defining the window in which I will look for negative vel for stim on turns
vel_xlim_rails_1=501; %%% defining the window in which I will look for negative vel for stim on rails
vel_xlim_rails_2=vel_xlim_rails_1+current_rails_dur; %%% defining the window in which I will look for negative vel for stim on rails

time_window_before =10*parameters.SampleRate;
time_window_after = (rails_dur_itr+additional_seconds_for_visualization)*parameters.SampleRate;
total_window_for_visualization=(2*time_window_after)+1;
total_time_for_visualization=rails_dur_itr+additional_seconds_for_visualization;
total_window_frames = time_window_before+time_window_after+1;

stimulus_intensities = [];

all_behavior_transitions_for_frame = {};
all_behavior_annotations_for_frame = {};
all_velocity_annotations_for_frame = {};
all_ellipse_ratio_annotations_for_frame = {};

%%%%% loop through each folder
velocity_at_stim_while_turning_folder=[];
velocity_at_stim_rails_folder=[];
ellipse_ratio_at_stim_while_turning_folder=[];
stim_while_turning_peaks_array_folder=[];
stim_on_rails_peaks_array_folder=[];
count_number_of_reversal_folder_turns=[];
count_number_of_reversal_folder_rails=[];
beh_state_on_stim_folder_new_protocol=[];  %%%% to find the turns on stim onset
ellipse_ratio_during_turn_state_folder=[];
just_the_ellipse_ratio_worm_data_folder=[];
possible_head_tail_switch_folder=[];
all_ellipse_ratio_folder=[];
cumulative_recording_duration_folder=0; %% cumulative reconding length
trial_count_folder=0;       %% to count the number of worm tracks in an assay
worms_tracked_simultaneously_folders=[];   %% to count the number of worms tracked simultaneously
error_code_folder=[]; %% to determine the criteria due to which a worm is not considered for analysis
width_vel_folder=[];
test=[];
fraction_occupancy_folder_turns=[];
fraction_occupancy_folder_rails=[];
behavior_on_additional_window_folder=[];
velocity_on_additional_window_folder=[];
ellipse_ratio_on_additional_window_folder=[];

for folder_index = [1:length(folders)]

% % % %     specific_folders_of_interest = [1:3];
% % % %     if ~ismember(folder_index,specific_folders_of_interest)
% % % %         continue
% % % %     end

    %%%% let us load the parameter file for each folder to determine the worm length threshold
    parameters_length = load_parameters(folders{folder_index});
    length_of_worm_threshold=minimum_length_of_worm*parameters_length.CameraPixeltommConversion;
    
    sprintf('Analyzing data %d out of %d datasets for %d sec duration', folder_index, length(folders),current_rails_dur/parameters.SampleRate)
    %%%%load the tracks for this folder
    [current_tracks, ~, ~] = loadtracks(folders{folder_index},relevant_track_fields);
     
    %%%%% generate the Behavior matricies
    current_tracks = BehavioralTransitionToBehavioralAnnotation(current_tracks);
    current_tracks = get_behavior_triggers(current_tracks);
    
    if size(current_tracks,2)==0  %%%% skip the iteration if it is an empty worm track
        continue
    end
        
    all_behavior_data{:,folder_index} = horzcat(current_tracks.BehavioralAnnotation);
    
    current_param = load_parameters(folders{folder_index});
    
    stim_while_turning_peaks_array_trial=[];
    stim_on_rails_peaks_array_trial=[];
    fullrails_peaks_array=[];
    
    velocity_at_stim_while_turning_trial=[];
    ellipse_ratio_at_stim_while_turning_trial=[];
    count_number_of_reversal_trial_turns=[];
    
    velocity_at_stim_while_turning_trial_rails=[];
    ellipse_ratio_at_stim_on_rails_trial=[];
    count_number_of_reversal_trial_rails=[];
    just_the_ellipse_ratio_worm_data_trial=[];
    possible_head_tail_switch_trial=[];

    beh_state_on_stim_trial_new_protocol=[];  %%%% to find the turns on stim onset
    ellipse_ratio_during_turn_state_trial=[];
    all_ellipse_ratio_trial=[];
    cumulative_recording_duration_trial=0;
    trial_count_trial=0;
    frame_details_trial=[];
    error_code_trial=[];
    pks_euc_dis_trial=[];
    pks_mean_angle_trial=[];
    pks_curvature_trial=[];
    width_vel_trial=[];
    all_valid_stim_rails_array_trial=[];
    fraction_occupancy_trial_turns=[];
    fraction_occupancy_trial_rails=[];
    behavior_on_additional_window_trial=[];
    velocity_on_additional_window_trial=[];
    ellipse_ratio_on_additional_window_trial=[];
    
    for track_index = 1:length(current_tracks)
%     for track_index=[183]
%     for track_index=[30 183 236 266 268 327 341 493 505 514 517 552]
%     for track_index = [60 91 183 236 292 327 341 414 415 416 420 426 484] %%%% [60 91 183 236 283 292]  %% for debugging
%     for track_index=590
        %%%% calculating cumulative recording duration and num of trials for each trial
        cumulative_recording_duration_trial=cumulative_recording_duration_trial+size(current_tracks(track_index).Frames,1);
        trial_count_trial=size(current_tracks,2);
        
        %%%% cancatenating Frame informaion to find out the number of worms tracked simultaneously at a time.
        
        frame_details_trial=[frame_details_trial
            current_tracks(track_index).Frames];
                
        %%%%% if no stimulus was delivered then ignore the iteration
        if nnz(current_tracks(track_index).AlignedStimulus(:,10))==0
            continue
        end
        
        if two_stim_delivered==1
           current_tracks(track_index).AlignedStimulus(:,:,1)=current_tracks(track_index).AlignedStimulus(:,:,opto_stim_color); 
        end
            
        mid_cline_index=size(current_tracks(track_index).AlignedStimulus,2)/2;
        
% % % %         %%%% if the stim is blue then multiply stim by -1closestIndex
% % % %         if min(min(current_tracks(track_index).AlignedStimulus))<0
% % % %             current_tracks(track_index).AlignedStimulus(:,:)=-1*current_tracks(track_index).AlignedStimulus(:,:); 
% % % %         end

        %%%% Whenever the stimulus is negative, make it positive
        for pp=1:size(current_tracks(track_index).AlignedStimulus,1)
            if current_tracks(track_index).AlignedStimulus(pp,:)<0
                current_tracks(track_index).AlignedStimulus(pp,:)=-1*current_tracks(track_index).AlignedStimulus(pp,:); 
            end
        end
                
        %%%%%%%%%%%%%%%%% to ignore smaller length worms %%%%%%%%%%%%%%%%%%
        length_of_current_worm=mean(current_tracks(track_index).Length);
        length_of_worm_matrix{folder_index,track_index}=length_of_current_worm;
        
        if length_of_current_worm<length_of_worm_threshold   %%% ignoring the lower 10 percentile worm length
            dummy_error_code(:,1)=folder_index;
            dummy_error_code(:,2)=track_index;
            dummy_error_code(:,3)=1;
            dummy_error_code(:,4)=1;
            error_code_trial=[error_code_trial
                dummy_error_code];    %%% error code for short length =1
            continue
        end

        %%%%%%%%%%%%%%% determining missing behavior indices %%%%%%%%%%%%%%
        
        original_velocity_for_correction=current_tracks(track_index).Velocity;
        smooth_velocity_for_correction=smooth(original_velocity_for_correction,10);   %%%% smoothing the velocity data
        current_tracks(track_index).Velocity=smooth_velocity_for_correction;  %%% redefinig the velocity info with smootherd data
        smooth_velocity_for_correction=smooth_velocity_for_correction';
        original_velocity_matrix_for_correction{track_index,:}=smooth_velocity_for_correction;
        
        original_ellipse_ratio_for_correction=current_tracks(track_index).EllipseRatio;
        smooth_ellipse_ratio_for_correction=smooth(original_ellipse_ratio_for_correction,10)'; %% smotting ellipse ratio
        current_tracks(track_index).EllipseRatio=smooth_ellipse_ratio_for_correction'; %% redefining the ellipse ratio with smoothed data
        smooth_ellipse_ratio_for_correction_matrix{track_index,:}=smooth_ellipse_ratio_for_correction;
        
        %%%%%% Finding the ellipse ratio when the worm is in the turn state
        index_during_turn_state=find(current_tracks(track_index).BehavioralAnnotation==2);
        ellipse_ratio_during_turn_state=current_tracks(track_index).EllipseRatio([index_during_turn_state],1);
        ellipse_ratio_during_turn_state_trial=[ellipse_ratio_during_turn_state_trial
            ellipse_ratio_during_turn_state];
        
        %%%%% cancatenating all the ellipse ratio of the worm
        all_ellipse_ratio_track=current_tracks(track_index).EllipseRatio;
        all_ellipse_ratio_trial=[all_ellipse_ratio_trial 
            all_ellipse_ratio_track];

        %%%% if pipeline does'nt assign behavior to one of the three states then we will take velocity into consideration to assign forward or reverse
        %%%% assigning beh as turn when the ellipse ratio was below the ellipse ration threshold
        for xyz=1:size(current_tracks(track_index).EllipseRatio,1)
            if smooth_ellipse_ratio_for_correction(1,xyz)<ellipse_ratio_threshold
                current_tracks(track_index).BehavioralAnnotation(1,xyz)=2;
            end
        end
        
        original_data_pre_correction=current_tracks(track_index).BehavioralAnnotation;
        original_data_matrix_pre_correction{track_index,:}=original_data_pre_correction;
        
        data_unassigned_pre_correction=find(original_data_pre_correction==0); %%%% when beh data is unassigned
        
        corrected_data_for_missing_behaviors=original_data_pre_correction;
        for xx=[data_unassigned_pre_correction]
            if smooth_velocity_for_correction(1,xx)<=0   %%%% -0.06
                corrected_data_for_missing_behaviors(1,xx)=3;   %%% assigning reversal
            else
                corrected_data_for_missing_behaviors(1,xx)=1;   %%% assigning forward
            end
        end
        
        current_tracks(track_index).BehavioralAnnotation=corrected_data_for_missing_behaviors;  %%%% we are reassigning behavior index for missing frames
        corrected_data_matrix_for_missing_behaviors{track_index,:}=corrected_data_for_missing_behaviors;
        
        %%% this is answer reviewer question on how reliable was the 10 s duration to detect head
        %%% I am finding the amount of time spent in reversal
        
        [~,~,width_vel,~]=findpeaks(double(smooth_velocity_for_correction<-0.06));
        dummy_width_vel=[];
        for xyz=1:size(width_vel,2)
            dummy_width_vel(xyz,1)=track_index;
            dummy_width_vel(xyz,2)=width_vel(1,xyz);
        end
        width_vel_trial=[width_vel_trial; dummy_width_vel];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        original_stim=current_tracks(track_index).AlignedStimulus(:,mid_cline_index);
        [raw_stim_peaks_amplitude, raw_stim_peaks_locs,raw_stim_widths,~] = findpeaks(current_tracks(track_index).AlignedStimulus(:,mid_cline_index), 'MinPeakDistance',current_param.InterTriggerInterval);
        
        %%%%BEGIN ANDYS REWRITE
        
        %%%start with the original stimulus time series and binarize the
        %%%stimulus into on and off depending on whether it exceeds a stimulu
        %%%sthreshold; then take the derivatives and find instances above zero 
        %%%and add 1 to account for teh derivative time shift effect
        stim_onsets = find(diff([original_stim>stim_threshold_min & original_stim<stim_threshold_max]) >0)+1;
        
        %%%Now we will apply increasingly selective criteria to only end up 
        %%%with stim_onset events that occur exactly during the onset of a turn
        
        %%%require that thes occur during times that mochi's pipeline classifies as turn
        data_turn=find(current_tracks(track_index).BehavioralAnnotation==behavior_to_study); %% when worm is turning
        stims_delivered_during_mochis_turns = intersect(stim_onsets, data_turn);
        
        %%%get the binary timeseries for when ellipse ratio dips below thresohold then take diff        
        er_neg_crossing_timeseries = [0; diff((current_tracks(track_index).EllipseRatio<ellipse_ratio_threshold))>0];
        
        %%%%then find negative values (which are the negative crossing events)
        er_turn_onset_indices = find(er_neg_crossing_timeseries>0);
        
        %%%Then get a list of all indices within +-10 frames of the crossing
        %%%event
        er_turn_onset_zones_indices = find(movmax((er_neg_crossing_timeseries>0),round(parameters.SampleRate/1.5))==1);
        
        %%%then filter valid stims further and impose that they fall in this
        %%%temporal window tied to the ellipsoid ratio crossing event
% % % %         valid_stim = intersect(stims_delivered_during_mochis_turns,er_turn_onset_zones_indices);
        valid_stim = stims_delivered_during_mochis_turns;
        valid_stim_matrix{track_index,:}=valid_stim;
        
% % % %         final_valid_stim=find_stim_on_turn_onsets(stims_delivered_during_mochis_turns,er_turn_onset_indices);
        final_valid_stim=stims_delivered_during_mochis_turns;
        final_valid_stim_matrix{track_index,:}=final_valid_stim;
    
        for mn=1:size(stim_onsets,1)
            
            beh_state_on_stim_track_new_protocol(:,1)=folder_index;
            beh_state_on_stim_track_new_protocol(:,2)=track_index;
            beh_state_on_stim_track_new_protocol(:,3)=stim_onsets(mn);
            
            if ismember(stim_onsets(mn),final_valid_stim)
                beh_state_on_stim_track_new_protocol(:,4)=1;  %% these are the valid stims on turn onset
            else
                beh_state_on_stim_track_new_protocol(:,4)=0;  %% these are the in-valid stims. They were not delivered on turn onset
            end
            
            beh_state_on_stim_trial_new_protocol=[beh_state_on_stim_trial_new_protocol
                beh_state_on_stim_track_new_protocol];    
            
        end

        %%% END ANDYS REWRITE
        
        raw_stim_peaks_amplitude=raw_stim_peaks_amplitude';
        raw_stim_peaks_locs=raw_stim_peaks_locs';
        raw_stim_widths=raw_stim_widths';
        
        %%%% determining peaks of the same widths                      
        raw_stim_peaks=[];
        raw_stim_peaks_height=[];
        raw_stim_duration=[];
        for i=1:size(raw_stim_widths,2)
            if (abs(current_rails_dur - raw_stim_widths(1,i)))<max_allowed_stim_duration_diff
            raw_stim_peaks=[raw_stim_peaks raw_stim_peaks_locs(1,i)];  %%% these are all the peaks with same width
            raw_stim_peaks_height=[raw_stim_peaks_height raw_stim_peaks_amplitude(1,i)];
            raw_stim_duration=[raw_stim_duration raw_stim_widths(1,i)];
            end
        end
        
        %%%%% consider stim of particular duration for this iteration
        
        new_stim_data=zeros(size(original_stim));
        for i=1:size(raw_stim_peaks,2)
            new_stim_data([raw_stim_peaks(1,i):(raw_stim_peaks(1,i)+raw_stim_duration(1,i))],1)=raw_stim_peaks_height(1,i);
        end
        
        new_stim_data(new_stim_data<stim_threshold_min)=0;
        new_stim_data(new_stim_data>stim_threshold_max)=0;
        
        if analysis_stim_while_turning==1  
        %%%% determining when the stim in on/off turn
        data_stim=find(new_stim_data>stim_threshold_min & new_stim_data<stim_threshold_max); %%%% when worm is stimulated

        raw_stim_final=zeros(1,size(current_tracks(track_index).AlignedStimulus,1));
        
        raw_data_stim_peaks=[];
        raw_stim_final(data_stim)=max(current_tracks(track_index).AlignedStimulus(:,mid_cline_index));
        
        [~, raw_data_stim_peaks] = findpeaks(raw_stim_final, 'MinPeakDistance',current_param.InterTriggerInterval); 
        
        all_data_stim_peaks_cell{track_index,1}=raw_data_stim_peaks;  %% we will capture all the raw stim peaks
        num_all_data_stim_peaks_array(track_index,1)=size(raw_data_stim_peaks,2);

        stim_on_turn=intersect(data_stim,data_turn); %%% correct stim based on turn 

        stim_on_turn_final=zeros(1,size(current_tracks(track_index).AlignedStimulus,1));
        stim_on_turn_final(stim_on_turn)=max(current_tracks(track_index).AlignedStimulus(:,mid_cline_index)); %%% this is the correct stimulation array
        [~, stim_on_turn_peaks] = findpeaks(stim_on_turn_final, 'MinPeakHeight',0.1, 'MinPeakDistance',current_param.InterTriggerInterval); 
        
        [minValue, closestIndex] = min(abs(stim_on_turn_peaks - raw_data_stim_peaks.'));
        stim_while_turning_peaks_dummy=unique(raw_data_stim_peaks(closestIndex));   %%%% unique is making sure that peaks do not get counted multiple times
        
        if behavior_to_study==2 || behavior_to_study==3
        stim_while_turning_peaks_old=[];
        
        for yz=1:length(stim_while_turning_peaks_dummy)
            if current_tracks(track_index).BehavioralAnnotation(1,stim_while_turning_peaks_dummy(yz))==behavior_to_study
                stim_while_turning_peaks_old=[stim_while_turning_peaks_old stim_while_turning_peaks_dummy(yz)];
            end
        end
        end
        
        stim_while_not_turning_peaks=setdiff(raw_data_stim_peaks,stim_while_turning_peaks_old);
        
        original_stim=current_tracks(track_index).AlignedStimulus(:,mid_cline_index);
        stim_while_turning_final=zeros(size(original_stim));
        
        stim_while_turning_peaks=intersect(final_valid_stim,raw_stim_peaks); %%% this array has the correctly defined stim associated turn peaks and peaks of same width 
        
        final_array_of_stim_while_turning_peaks=[];    %%%% collecting just the final stim instances which goes through the analysis pipeline. It is for visualization code
        for i=1:length(stim_while_turning_peaks)
            
            %%%% sprintf('current_track_index: %d and stim_at_peaks: %d',track_index,stim_while_turning_peaks(i))
            
            %%% Noting down all the stim while turning peaks array
            dummy_stim_on_valid_turns(:,1)=folder_index;
            dummy_stim_on_valid_turns(:,2)=track_index;
            dummy_stim_on_valid_turns(:,3)=stim_while_turning_peaks(i);
            
            stim_while_turning_peaks_array_trial=[stim_while_turning_peaks_array_trial
                dummy_stim_on_valid_turns];
                        
            %%% to avoid error due to 'findpeaks' function when the array begins with peak
            if any(stim_while_turning_peaks(i)<10)
                %%% disp('Error 10: Stim on turn peaks within 10 frames')
                dummy_error_code(:,1)=folder_index;
                dummy_error_code(:,2)=track_index;
                dummy_error_code(:,3)=[];
                dummy_error_code(:,4)=9; %% error code is 9 when worm track begins with stim
                error_code_trial=[error_code_trial
                    dummy_error_code];  
                continue
            end
        
            if stim_while_turning_peaks(i)<=500
                dummy_error_code(:,1)=folder_index;
                dummy_error_code(:,2)=track_index;
                dummy_error_code(:,3)=stim_while_turning_peaks(i);
                dummy_error_code(:,4)=2;
                error_code_trial=[error_code_trial
                dummy_error_code];    %%% if the worm is not tracked for 17 sec pre-stimulus onset error code =2
                continue
            end
            
            if size(current_tracks(track_index).AlignedStimulus(:,10),1)<=(stim_while_turning_peaks(i)+current_rails_dur+(additional_seconds_for_visualization*parameters.SampleRate))
                dummy_error_code(:,1)=folder_index;
                dummy_error_code(:,2)=track_index;
                dummy_error_code(:,3)=stim_while_turning_peaks(i);
                dummy_error_code(:,4)=3;
                error_code_trial=[error_code_trial
                dummy_error_code];    %%% error code=3 if the worm is not tracked for the entirety stimulus duration
                continue
            end
                
            dummy_array_with_velocity_info_turns=[];
            dummy_array_with_ellipse_ratio_info_turns=[];
            dummy_array_with_ellipse_ratio_info_turns=[];
            locs_er=[];
            time_spend_in_turns=[];
            
            index_turn_peak=find(stim_while_turning_peaks(i)==raw_stim_peaks_locs);
            stim_while_turning_final([raw_stim_peaks_locs(1,index_turn_peak):(raw_stim_peaks_locs(1,index_turn_peak)+raw_stim_widths(1,index_turn_peak))],1)=raw_stim_peaks_amplitude(1,index_turn_peak);
            
            dummy_array_with_velocity_info_turns=current_tracks(track_index).Velocity((stim_while_turning_peaks(i)-500):(stim_while_turning_peaks(i)+current_rails_dur),1);
            
            velocity_at_stim_while_turning_trial=[velocity_at_stim_while_turning_trial
                dummy_array_with_velocity_info_turns];
            
            if max(dummy_array_with_velocity_info_turns(1:400,:))<=0.01 && min(dummy_array_with_velocity_info_turns(1:400,:))>=-0.01  %%%% ignore trials when worms hardly moved
                %%% disp('Worm is not moving: Trial ignored')
                dummy_error_code(:,1)=folder_index;
                dummy_error_code(:,2)=track_index;
                dummy_error_code(:,3)=stim_while_turning_peaks(i);
                dummy_error_code(:,4)=4;
                error_code_trial=[error_code_trial
                dummy_error_code];    %%% error code for not moving worm =4
                continue
            end
            
            dummy_array_with_ellipse_ratio_info_turns=current_tracks(track_index).EllipseRatio((stim_while_turning_peaks(i)-500):(stim_while_turning_peaks(i)+current_rails_dur),1);
       
            just_the_ellipse_ratio_worm_data_trial=[just_the_ellipse_ratio_worm_data_trial
                dummy_array_with_ellipse_ratio_info_turns];
            
            %%%%%%%%%% to get rid of trials when two worms collide
            if any(dummy_array_with_ellipse_ratio_info_turns>8)
                %%% disp('worms collided during trial')
                dummy_error_code(:,1)=folder_index;
                dummy_error_code(:,2)=track_index;
                dummy_error_code(:,3)=stim_while_turning_peaks(i);
                dummy_error_code(:,4)=5;
                error_code_trial=[error_code_trial
                dummy_error_code];    %%% error code when worm collide =5
                continue
            end
            
            %%%%% to get rid of trials with very low ellipse ratio in turns analysis 
            rng(1,'twister')
            randIdcs_turns = randi([1 400],1,50);
            test_low_ellipse_ratio_worm = dummy_array_with_ellipse_ratio_info_turns(randIdcs_turns,:);

            if mean(test_low_ellipse_ratio_worm)<test_low_ellipse_ratio_worm_threshold
                %%% disp('Worm with low mean ellipse ratio')
                dummy_error_code(:,1)=folder_index;
                dummy_error_code(:,2)=track_index;
                dummy_error_code(:,3)=stim_while_turning_peaks(i);
                dummy_error_code(:,4)=6;
                error_code_trial=[error_code_trial
                dummy_error_code];    %%% error code when worm's mean ellipse ratio is low =6
                continue
            end
        
            ellipse_ratio_at_stim_while_turning_trial=[ellipse_ratio_at_stim_while_turning_trial
                dummy_array_with_ellipse_ratio_info_turns];
            
% % %             %%% make sure that stim is delivered during turns and not reversals 
% % %             if any(dummy_array_with_velocity_info_turns(495:500,:)<(negative_vel_threshold/2))
% % %                 %%% disp('Stim while reversal: trial ignored')
% % %                 dummy_error_code(:,1)=folder_index;
% % %                 dummy_error_code(:,2)=track_index;
% % %                 dummy_error_code(:,3)=stim_while_turning_peaks(i);
% % %                 dummy_error_code(:,4)=7;
% % %                 error_code_trial=[error_code_trial
% % %                 dummy_error_code];    %%% error code when stim delivered to reversing worms =7
% % %                 continue
% % %             end
            
            %%% detecting head-tail switch
            current_mean_angles=current_tracks(track_index).MeanAngle(:,(stim_while_turning_peaks(i)-1):(stim_while_turning_peaks(i)+current_rails_dur));
            centerlines_data=current_tracks(track_index).Centerlines(:,:,(stim_while_turning_peaks(i)-1:(stim_while_turning_peaks(i)+current_rails_dur)));
            
            [possible_head_tail_switch_track,pks_euc_dis_track,pks_mean_angle_track,pks_curvature_track]=detect_head_tail_switch_3(current_mean_angles, centerlines_data,stim_while_turning_peaks(i),current_rails_dur,folder_index,track_index,plot_head_tail_switch_instances);
            possible_head_tail_switch_trial=[possible_head_tail_switch_trial; possible_head_tail_switch_track];
            
            pks_euc_dis_trial=[pks_euc_dis_trial pks_euc_dis_track'];
            pks_mean_angle_trial=[pks_mean_angle_trial pks_mean_angle_track];
            pks_curvature_trial=[pks_curvature_trial pks_curvature_track];
            
            %%%%%% to detect if the worm keeps on reversing for a longer time pre-stimulus
            data_reverse_correction_turns=current_tracks(track_index).BehavioralAnnotation(:,stim_while_turning_peaks(i)-500:stim_while_turning_peaks(i)-100); %%%% when worm is reversing

            [~, reverse_locs_correction_turns,reverse_widths_correction_turns,~] = findpeaks(data_reverse_correction_turns,'MinPeakHeight',2.5);
            if any(reverse_widths_correction_turns>210) 
                %%% disp('Head tail switch pre stimulus')
                dummy_error_code(:,1)=folder_index;
                dummy_error_code(:,2)=track_index;
                dummy_error_code(:,3)=stim_while_turning_peaks(i);
                dummy_error_code(:,4)=10;
                error_code_trial=[error_code_trial
                dummy_error_code];    %%% error code when head tail switch occurs prestimuls =70
                continue 
            end
            
            %%%%% to generate a plot for fraction occpancy during turns
            fraction_occupancy_dummy_turns_labels(:,1)=folder_index;
            fraction_occupancy_dummy_turns_labels(:,2)=track_index;
            fraction_occupancy_dummy_turns_labels(:,3)=stim_while_turning_peaks(i);
            fraction_occupancy_dummy_turns=current_tracks(track_index).BehavioralAnnotation(1,stim_while_turning_peaks(i)-time_window_before:stim_while_turning_peaks(i)+time_window_after);   
            fraction_occupancy_dummy_turns_with_labels=[fraction_occupancy_dummy_turns_labels fraction_occupancy_dummy_turns];
            
            fraction_occupancy_trial_turns=[fraction_occupancy_trial_turns
                fraction_occupancy_dummy_turns_with_labels];
            
            %%%% collecting all the valid peaks for visualization
            final_array_of_stim_while_turning_peaks=[final_array_of_stim_while_turning_peaks
                stim_while_turning_peaks(i)];
            
            %%%%%% detecting reversals in the stimulus window 
            indices_below_set_vel_threshold_turns = find(dummy_array_with_velocity_info_turns(vel_xlim_turns_1:vel_xlim_turns_2,:) < negative_vel_threshold); 

            if size(indices_below_set_vel_threshold_turns,1)>=frames_in_reversal_threshold 
                dummy_array_with_number_of_reversal_turns(:,1)=folder_index;
                dummy_array_with_number_of_reversal_turns(:,2)=track_index;
                dummy_array_with_number_of_reversal_turns(:,3)=stim_while_turning_peaks(i);
                dummy_array_with_number_of_reversal_turns(:,4)=1;
            else
                dummy_array_with_number_of_reversal_turns(:,1)=folder_index;
                dummy_array_with_number_of_reversal_turns(:,2)=track_index;
                dummy_array_with_number_of_reversal_turns(:,3)=stim_while_turning_peaks(i);
                dummy_array_with_number_of_reversal_turns(:,4)=0;
            end

            count_number_of_reversal_trial_turns=[count_number_of_reversal_trial_turns
                dummy_array_with_number_of_reversal_turns];  %%% cancatenating the info if a worm reversed or not
                                               
        end
        
        if analysis_stim_while_turning==1
            single_dimension_stimulus=stim_while_turning_final';
        end
        
        end
        
        %%%%%%%%%%%%%%%%%%% analysis for open loop %%%%%%%%%%%%%%%%%%%%%%%%
        
        if analysis_rails==1
            single_dimension_stimulus=new_stim_data';
        end
        
        %%%%% all stim rails of valid duration and intensity
        all_valid_stim_rails_array_trial{track_index,1}=single_dimension_stimulus;
        
        xcorr_ledvoltages_stimulus = abs(padded_conv(single_dimension_stimulus, normalized_stimuli));
        [~, all_stim_peaks] = findpeaks(xcorr_ledvoltages_stimulus, 'MinPeakDistance',current_param.InterTriggerInterval);

        %%%%%%%%%%%% counting number of reversals for rails %%%%%%%%%%%%%%%
        final_array_of_stim_on_rails_peaks=[];
        if analysis_rails==1
            stim_peaks=setdiff(all_stim_peaks,final_valid_stim); %%% get rid of the stims which occured on the onset of turns
            
            for ab=1:length(stim_peaks)
            
            %%% sprintf('current_track_index: %d and stim_at_peaks: %d',track_index,stim_peaks(ab))
            
            %%% Noting down all the stim on rails
            dummy_stim_on_rails(:,1)=folder_index;
            dummy_stim_on_rails(:,2)=track_index;
            dummy_stim_on_rails(:,3)=stim_peaks(ab);
            
            stim_on_rails_peaks_array_trial=[stim_on_rails_peaks_array_trial
                dummy_stim_on_rails];
            
            if stim_peaks(ab)<=500
                %%% disp('Stim in the first 500 frames: trial ignored')
                dummy_error_code(:,1)=folder_index;
                dummy_error_code(:,2)=track_index;
                dummy_error_code(:,3)=stim_peaks(ab);
                dummy_error_code(:,4)=2;   %% worm not tracked for 17 seconds before stimulus onset
                error_code_trial=[error_code_trial 
                dummy_error_code];   
                continue
            end
            
            if size(current_tracks(track_index).AlignedStimulus(:,10),1)<=(stim_peaks(ab)+current_rails_dur+(additional_seconds_for_visualization*parameters.SampleRate))
                dummy_error_code(:,1)=folder_index;
                dummy_error_code(:,2)=track_index;
                dummy_error_code(:,3)=stim_peaks(ab);
                dummy_error_code(:,4)=3; %% worm not tracked for the complete stimulus duration
                error_code_trial=[error_code_trial
                dummy_error_code];  
                continue
            end
               
            dummy_array_with_velocity_info_rails=[];
            
            stim_on_rails_peaks_array(track_index,ab)=stim_peaks(ab);
            
% % % % %             %%%%%% this portion is only to answer reviewer comments.
% % % % %             if stim_peaks(ab)+450>size(current_tracks(track_index).Frames,1)
% % % % %                 continue
% % % % %             end
            %%% i am making an array of velocity for 15 seconds  post stim dur. Originally
            %%% i made an array only for stim dur
            dummy_array_with_velocity_info_rails=current_tracks(track_index).Velocity((stim_peaks(ab)-500):(stim_peaks(ab)+450),1);
            
%             dummy_array_with_velocity_info_rails=current_tracks(track_index).Velocity((stim_peaks(ab)-500):(stim_peaks(ab)+current_rails_dur),1);
            
            if max(dummy_array_with_velocity_info_rails(1:400,:))<=0.01 && min(dummy_array_with_velocity_info_rails(1:400,:))>=-0.01  %%%% ignore trials when worms hardly moved
                %%% disp('Worm is not moving: Trial ignored')
                dummy_error_code(:,1)=folder_index;
                dummy_error_code(:,2)=track_index;
                dummy_error_code(:,3)=stim_peaks(ab);
                dummy_error_code(:,4)=4;  %% worm is not moving (stationary most of the time)
                error_code_trial=[error_code_trial
                dummy_error_code];  
                continue
            end
            
            velocity_at_stim_while_turning_trial_rails=[velocity_at_stim_while_turning_trial_rails
                dummy_array_with_velocity_info_rails'];
            
            dummy_array_with_ellipse_ratio_info_rails=current_tracks(track_index).EllipseRatio((stim_peaks(ab)-500):(stim_peaks(ab)+current_rails_dur),1);
            
            %%%%%%%%%% to get rid of trials when two worms collide
            if any(dummy_array_with_ellipse_ratio_info_rails>8)
                %%% disp('worms collided during trial')
                dummy_error_code(:,1)=folder_index;
                dummy_error_code(:,2)=track_index;
                dummy_error_code(:,3)=stim_peaks(ab);
                dummy_error_code(:,4)=5;
                error_code_trial=[error_code_trial
                dummy_error_code];    %%% error code when worm collide =5
                continue
            end
            
            %%%%% to get rid of trials with very low ellipse ratio in turns analysis 
            rng(1,'twister')
            randIdcs_rails = randi([1 400],1,50);
            test_low_ellipse_ratio_worm_rails = dummy_array_with_ellipse_ratio_info_rails(randIdcs_rails,:);
                     
            if mean(test_low_ellipse_ratio_worm_rails)<test_low_ellipse_ratio_worm_threshold
                %%% disp('Worm with low mean ellipse ratio')
                dummy_error_code(:,1)=folder_index;
                dummy_error_code(:,2)=track_index;
                dummy_error_code(:,3)=stim_peaks(ab);
                dummy_error_code(:,4)=6;  %% low mean ellipse ratio of worm
                error_code_trial=[error_code_trial
                dummy_error_code];  
                continue
            end

            ellipse_ratio_at_stim_on_rails_trial=[ellipse_ratio_at_stim_on_rails_trial
                dummy_array_with_ellipse_ratio_info_rails];
            
            %%% make sure that stim is delivereed during forward motion and not reversals
            if any(dummy_array_with_velocity_info_rails(495:500,:)<negative_vel_threshold/2)
                %%% disp('Stim while reversal: trial ignored')
                dummy_error_code(:,1)=folder_index;
                dummy_error_code(:,2)=track_index;
                dummy_error_code(:,3)=stim_peaks(ab);
                dummy_error_code(:,4)=7;  %% stim delivered during reversal
                error_code_trial=[error_code_trial
                dummy_error_code];  
                continue
            end
            
            %%% make sure that stim is delivered during forward motion and not turns
            if any(dummy_array_with_ellipse_ratio_info_rails(490:500,:)<=ellipse_ratio_threshold)
                %%% disp('Stim on turn: trial ignored')
                dummy_error_code(:,1)=folder_index;
                dummy_error_code(:,2)=track_index;
                dummy_error_code(:,3)=stim_peaks(ab);
                dummy_error_code(:,4)=8; %% stim delivered during turn
                error_code_trial=[error_code_trial
                dummy_error_code];  
                continue
            end
            
            %%%%%% to detect if the worm keeps on reversing for a longer time pre-stimulus
            data_reverse_correction_rails=current_tracks(track_index).BehavioralAnnotation(:,stim_peaks(ab)-500:stim_peaks(ab)-100); %%%% when worm is reversing

            [~, reverse_locs_correction_rails,reverse_widths_correction_rails,~] = findpeaks(data_reverse_correction_rails);

            if any(reverse_widths_correction_rails>210) 
                %%% disp('Head tail switch pre stimulus')
                dummy_error_code(:,1)=folder_index;
                dummy_error_code(:,2)=track_index;
                dummy_error_code(:,3)=stim_peaks(ab);
                dummy_error_code(:,4)=10;
                error_code_trial=[error_code_trial
                dummy_error_code];    %%% error code when head tail switch occurs prestimuls =70
                continue 
            end
            
            %%%%% to generate a plot for fraction occpancy during turns
            fraction_occupancy_dummy_rails=current_tracks(track_index).BehavioralAnnotation(1,stim_peaks(ab)-time_window_before:stim_peaks(ab)+time_window_after);   
            fraction_occupancy_trial_rails=[fraction_occupancy_trial_rails
                fraction_occupancy_dummy_rails];
            
            %%%%% collecting all the stim peaks for visualization plots
            final_array_of_stim_on_rails_peaks=[final_array_of_stim_on_rails_peaks
                stim_peaks(ab)];
            
           %%% here I am adding NaNs in case some tracks are not more than 3 sec longer since stim initiation 
           size_of_current_behavior_track=size(current_tracks(track_index).BehavioralAnnotation,2); 
           if (size_of_current_behavior_track-stim_peaks(ab))<=time_window_after 
               
               %%%% to debug
               test=[test; track_index];
                   
               original_size=size(current_tracks(track_index).BehavioralAnnotation,2);
               nans_to_add=time_window_after-(original_size-stim_peaks(ab));
               behavior_on_additional_window_track=[current_tracks(track_index).BehavioralAnnotation(1,stim_peaks(ab)-time_window_before:end) NaN(1,nans_to_add)];
               velocity_on_additional_window_track=[smooth_velocity_for_correction(1,stim_peaks(ab)-time_window_before:end) NaN(1,nans_to_add)];
               ellipse_ratio_on_additional_window_track=[smooth_ellipse_ratio_for_correction(1,stim_peaks(ab)-time_window_before:end) NaN(1,nans_to_add)];
           else
               behavior_on_additional_window_track=current_tracks(track_index).BehavioralAnnotation(1,stim_peaks(ab)-time_window_before:stim_peaks(ab)+time_window_after);
               velocity_on_additional_window_track=smooth_velocity_for_correction(1,stim_peaks(ab)-time_window_before:stim_peaks(ab)+time_window_after);
               ellipse_ratio_on_additional_window_track=smooth_ellipse_ratio_for_correction(1,stim_peaks(ab)-time_window_before:stim_peaks(ab)+time_window_after);
           end
           
          behavior_on_additional_window_trial=[behavior_on_additional_window_trial
               behavior_on_additional_window_track];
          velocity_on_additional_window_trial=[velocity_on_additional_window_trial
               velocity_on_additional_window_track];
          ellipse_ratio_on_additional_window_trial=[ellipse_ratio_on_additional_window_trial
               ellipse_ratio_on_additional_window_track];
           
            %%%% detecting reversals
            indices_below_set_vel_threshold_rails = find(dummy_array_with_velocity_info_rails(vel_xlim_rails_1:vel_xlim_rails_2,:) < negative_vel_threshold); 

            %%% when no reversal is detected
            if isempty(indices_below_set_vel_threshold_rails)

                %%% disp('No reversal detected')
                dummy_array_with_number_of_reversal_rails(:,1)=folder_index;
                dummy_array_with_number_of_reversal_rails(:,2)=track_index;
                dummy_array_with_number_of_reversal_rails(:,3)=stim_peaks(ab);
                dummy_array_with_number_of_reversal_rails(:,4)=0;
                dummy_array_with_number_of_reversal_rails(:,5)=current_tracks(track_index).Frames(stim_peaks(ab));
                
                count_number_of_reversal_trial_rails=[count_number_of_reversal_trial_rails
                dummy_array_with_number_of_reversal_rails];
                continue
            end
            
            %%%% if reversal is detected

            if size(indices_below_set_vel_threshold_rails,1)>=frames_in_reversal_threshold 
                dummy_array_with_number_of_reversal_rails(:,1)=folder_index;
                dummy_array_with_number_of_reversal_rails(:,2)=track_index;
                dummy_array_with_number_of_reversal_rails(:,3)=stim_peaks(ab);
                dummy_array_with_number_of_reversal_rails(:,4)=1;
                dummy_array_with_number_of_reversal_rails(:,5)=current_tracks(track_index).Frames(stim_peaks(ab));
            else
                dummy_array_with_number_of_reversal_rails(:,1)=folder_index;
                dummy_array_with_number_of_reversal_rails(:,2)=track_index;
                dummy_array_with_number_of_reversal_rails(:,3)=stim_peaks(ab);
                dummy_array_with_number_of_reversal_rails(:,4)=0;
                dummy_array_with_number_of_reversal_rails(:,5)=current_tracks(track_index).Frames(stim_peaks(ab));
            end

            count_number_of_reversal_trial_rails=[count_number_of_reversal_trial_rails
                dummy_array_with_number_of_reversal_rails];
                
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%% Mochi's code for visualization
        if analysis_rails==1
            stim_peaks=final_array_of_stim_on_rails_peaks;
        else
            stim_peaks=final_array_of_stim_while_turning_peaks;
        end
            
% % % % %         stim_peaks=all_stim_peaks;
        for peak_index = 1:length(stim_peaks)
            
            current_stim_power = (single_dimension_stimulus(stim_peaks(peak_index))); %%% i am rounding the stim power so that it does not accumulate zeros
            %%%%%% see if this stim_power already exists
            current_stim_index = find(or(and(stimulus_intensities <= current_stim_power*stim_similarity_thresh, stimulus_intensities >= current_stim_power/stim_similarity_thresh), ...
                and(-stimulus_intensities <= current_stim_power*-stim_similarity_thresh, -stimulus_intensities >= current_stim_power/-stim_similarity_thresh)));
            
            if isempty(current_stim_index)
                %%%%%% no entry yet, add it
                if any(raw_stim_peaks_height<=2) %%%%%% any(raw_stim_peaks_height<=2)
                stimulus_intensities = [stimulus_intensities,(current_stim_power)]; %%%%% when power is less than 1 then it accumulates zeros
                else
                stimulus_intensities = [stimulus_intensities,round(current_stim_power)];
                end
                
                current_stim_index = length(stimulus_intensities);
                all_behavior_transitions_for_frame{current_stim_index} = cell(1,total_window_frames);
                all_behavior_annotations_for_frame{current_stim_index} = cell(1,total_window_frames);
                all_velocity_annotations_for_frame{current_stim_index} = cell(1,total_window_frames);
                all_ellipse_ratio_annotations_for_frame{current_stim_index} = cell(1,total_window_frames);
            end
                
            for frame_shift = -time_window_before:time_window_after
                current_frame = stim_peaks(peak_index) + frame_shift;
                if current_frame <= length(current_tracks(track_index).Frames) && current_frame >= 1
                    %%%%%% make sure the current frame is in range
                    %%%%%% cut up tracks to each frame
                    all_behavior_transitions_for_frame{current_stim_index}{frame_shift+time_window_before+1} = [all_behavior_transitions_for_frame{current_stim_index}{frame_shift+time_window_before+1}, current_tracks(track_index).Behaviors(:,current_frame)];
                    all_behavior_annotations_for_frame{current_stim_index}{frame_shift+time_window_before+1} = [all_behavior_annotations_for_frame{current_stim_index}{frame_shift+time_window_before+1}, current_tracks(track_index).BehavioralAnnotation(current_frame)];
                    all_velocity_annotations_for_frame{current_stim_index}{frame_shift+time_window_before+1} = [all_velocity_annotations_for_frame{current_stim_index}{frame_shift+time_window_before+1}, current_tracks(track_index).Velocity(current_frame)];
                    all_ellipse_ratio_annotations_for_frame{current_stim_index}{frame_shift+time_window_before+1} = [all_ellipse_ratio_annotations_for_frame{current_stim_index}{frame_shift+time_window_before+1}, current_tracks(track_index).EllipseRatio(current_frame)];
                end
            end
        end
    end  
    
    %%%% counting the number of occurence of unique elements in frame_details_trial
    %%%% This will give me the number of worms tracked at each given frame
    [worms_tracked_simultaneously_trials,GroupR_frames_trials] = groupcounts(frame_details_trial);
    
    worms_tracked_simultaneously_folders=[worms_tracked_simultaneously_folders
        worms_tracked_simultaneously_trials];
    
    %%%% cancatenating the data from each assay (or trial) and generate for all the folders
    just_the_ellipse_ratio_worm_data_folder=[just_the_ellipse_ratio_worm_data_folder
        just_the_ellipse_ratio_worm_data_trial];
        
    velocity_at_stim_while_turning_folder=[velocity_at_stim_while_turning_folder
        velocity_at_stim_while_turning_trial];
    
    velocity_at_stim_rails_folder=[velocity_at_stim_rails_folder
    velocity_at_stim_while_turning_trial_rails];
    
    ellipse_ratio_at_stim_while_turning_folder=[ellipse_ratio_at_stim_while_turning_folder
    ellipse_ratio_at_stim_while_turning_trial];
    
    count_number_of_reversal_folder_turns=[count_number_of_reversal_folder_turns
        count_number_of_reversal_trial_turns];
    
    test=[test; size(count_number_of_reversal_trial_turns,1)];
    
    count_number_of_reversal_folder_rails=[count_number_of_reversal_folder_rails
    count_number_of_reversal_trial_rails];

    stim_while_turning_peaks_array_folder=[stim_while_turning_peaks_array_folder
        stim_while_turning_peaks_array_trial];
    
    stim_on_rails_peaks_array_folder=[stim_on_rails_peaks_array_folder
    stim_on_rails_peaks_array_trial];
    
    possible_head_tail_switch_folder=[possible_head_tail_switch_folder
        possible_head_tail_switch_trial];
    
    beh_state_on_stim_folder_new_protocol=[beh_state_on_stim_folder_new_protocol
        beh_state_on_stim_trial_new_protocol];
    
    ellipse_ratio_during_turn_state_folder=[ellipse_ratio_during_turn_state_folder
        ellipse_ratio_during_turn_state_trial];
    
    all_ellipse_ratio_folder=[all_ellipse_ratio_folder 
        all_ellipse_ratio_trial];
    
    error_code_folder=[error_code_folder
                error_code_trial];  
            
    width_vel_folder=[width_vel_folder
        width_vel_trial];
    
    fraction_occupancy_folder_turns=[fraction_occupancy_folder_turns
        fraction_occupancy_trial_turns];
    
    fraction_occupancy_folder_rails=[fraction_occupancy_folder_rails
        fraction_occupancy_trial_rails];
    
      behavior_on_additional_window_folder=[behavior_on_additional_window_folder
           behavior_on_additional_window_trial];
      velocity_on_additional_window_folder=[velocity_on_additional_window_folder
           velocity_on_additional_window_trial];
      ellipse_ratio_on_additional_window_folder=[ellipse_ratio_on_additional_window_folder
           ellipse_ratio_on_additional_window_trial];
    
    cumulative_recording_duration_folder=cumulative_recording_duration_folder+cumulative_recording_duration_trial;
    trial_count_folder=trial_count_folder+trial_count_trial;

end

if plotting_graphs==1 %%% to prevent graphs and corresponding calculations
%%%%%% 
%%% sort the stimulus intensities
[stimulus_intensities, sort_index] = sort((stimulus_intensities));

all_behavior_transitions_for_frame = all_behavior_transitions_for_frame(sort_index);
all_behavior_annotations_for_frame= all_behavior_annotations_for_frame(sort_index);
all_velocity_annotations_for_frame= all_velocity_annotations_for_frame(sort_index);
all_ellipse_ratio_annotations_for_frame= all_ellipse_ratio_annotations_for_frame(sort_index);

stimulus_intensities=round(stimulus_intensities);
n_sti=length(stimulus_intensities);
behavior_counts_for_frame = zeros(number_of_behaviors,n_sti,total_window_frames);
behavior_ratios_for_frame = zeros(number_of_behaviors,n_sti,total_window_frames);

for stimulus_index = 1:n_sti
    %%%%%% plot the transition rates centered on stim delivery
    total_counts_for_frame = zeros(1,total_window_frames);
    for frame_index = 1:total_window_frames
        for behavior_index = 1:number_of_behaviors
            behavior_counts_for_frame(behavior_index,stimulus_index,frame_index) = sum(all_behavior_annotations_for_frame{stimulus_index}{frame_index}==behavior_index);
        end
        behavior_ratios_for_frame(:,stimulus_index,frame_index) = behavior_counts_for_frame(:,stimulus_index,frame_index)./sum(behavior_counts_for_frame(:,stimulus_index,frame_index)); %get ratio
    end
end

%% 2 plot the behavioral ratios as a function of time

n_tracks=zeros(1,n_sti); %number of tracks in each sti intensities
for stimulus_index = 1:n_sti
    n_tracks(stimulus_index) = round(mean(arrayfun(@(x) size(x{1},2), [all_behavior_transitions_for_frame{1,stimulus_index}])));
    my_colors = behavior_colors;
    
    if plot_all_beh_ratios==1
    figure
    light_pulse = area([0 current_rails_dur/parameters.SampleRate], [1 1],'facecolor',stim_color,'LineStyle','none');
    alpha(0.1)
    hold on;
    for behavior_index = 1:number_of_behaviors
        plot(-time_window_before/fps:1/fps:time_window_after/fps, squeeze(behavior_ratios_for_frame(behavior_index,stimulus_index,:)), '-', 'color', my_colors(behavior_index,:),'Linewidth', 3,'DisplayName',behavior_names{behavior_index});
    end
    hold off
    xlabel('Time (s)') % x-axis label
    ylabel('Behavioral Ratio') % y-axis label
    ax = gca;
    
    title(['Int. = ', num2str(stimulus_intensities(stimulus_index)),'\muW/mm^2' ...
        ',' ' Width = ', num2str(current_rails_dur/parameters.SampleRate),'s' ...
        ',' ' n = ', num2str(n_tracks(stimulus_index)), ' events']);
    ax.FontSize = 13;
    %     xline(-6,'b--','DisplayName','LED on');xline(2,'b-.','DisplayName','LED off');
    legend('show','Location','northwest','FontSize',12);
%     tap_line=xline(0,'k--');
    set(get(get(light_pulse,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
%     set(get(get(tap_line,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    axis([-time_window_before/fps time_window_after/fps 0 1])
    box off;
%     ylim([0 0.3])
    end
%     axis([-time_window_before/fps time_window_after/fps 0 1])
end

%% plot the behavioral ratios for various intensities on the same plot
if plot_various_int==1
    
for behavior_index = [1,2,3] %fast forward3; fast reverse; turns
    my_colors = lines(n_sti);
    
    figure
    light_pulse = area([0 current_rails_dur/parameters.SampleRate], [1 1],'facecolor',stim_color,'LineStyle','none');
    alpha(0.1)
    hold on; 
    for stimulus_index = 1:n_sti
        %track_n = round(mean(arrayfun(@(x) size(x{1},2), [all_behavior_transitions_for_frame{1,stimulus_index}])));
        plot(-time_window_before/fps:1/fps:time_window_after/fps, squeeze(behavior_ratios_for_frame(behavior_index,stimulus_index,:)), '-', 'color', my_colors(stimulus_index,:),'Linewidth', 3,...
            'DisplayName',[num2str(stimulus_intensities(stimulus_index)), '\muW (n = ', num2str(n_tracks(stimulus_index)),')']);
    individual_behavior_data{behavior_index,stimulus_index}=squeeze(behavior_ratios_for_frame(behavior_index,stimulus_index,:));
    end
    hold off
    xlabel('Time (s)') % x-axis label
    ylabel([behavior_names{behavior_index}, ' Ratio'])
    title(['Beh. = ', behavior_names{behavior_index},...
        ',' ' Stim width = ', num2str(current_rails_dur/parameters.SampleRate),' s']);
    
%     xline(-2,'b--','DisplayName','LED on');xline(2,'b-.','DisplayName','LED off');
    set(get(get(light_pulse,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    ax = gca;
    ax.FontSize = 14;
    legend('show','Location','northwest','FontSize',12);
%     ylim([0 0.3])
%     yticks([0:0.1:0.3])
    xlim([-time_window_before/fps time_window_after/fps])
    xticks([-10:5:12])
    alpha(0.2)
    box off;
    grid off;
end
end
%% plotting mean and std error for a specific behavior state

if plot_comparison_reverse_beh==1
%optotap_behavioral_ratio_percent_changes = percent_change_above_baseline(squeeze(behavior_ratios_for_frame(:,4,:)));
behavior_index = [2,3];

% mn=behavior_index;
for mn=behavior_index

    % [behavior_percent_change, behavior_baselines, behavior_max, behavior_min] = percent_change_above_baseline(squeeze(behavior_ratios_for_frame(behavior_index,:,:)));
    % y = behavior_max';
    stimulus_intensity_array = stimulus_intensities;
    err_boot=zeros(1,length(stimulus_intensity_array));
    % standard error of the mean
    sem=zeros(1,length(stimulus_intensity_array));
    % find the max ratio within the time window
    [y2,mi]= max(behavior_ratios_for_frame(mn,:,301:end),[],3); %%%% To get values in specific range change the range of third element e.g. (behavior_index,:,301:end),[],3)
    for stimulus_index = 1:n_sti
        % create binary data at the frame when the ratio of behavior is max
% % % %         binary_data=(all_behavior_annotations_for_frame{stimulus_index}{mi(stimulus_index)}==mn);
        binary_data=(all_behavior_annotations_for_frame{stimulus_index}{mi(stimulus_index)+301}==mn);  %%%% Adding the same index number is very important
        binary_data_matrix{mn,:}=binary_data;
        % get standard error from bootstraped data
        % bootratio=bootstrp(20000,@mean,binary_data);
        % err_boot(stimulus_index)=std(bootratio);
        % mathematical formula of standard error of the mean
        sem(stimulus_index)=std(binary_data)/sqrt(length(binary_data));
    end
    
    if mn==2
        turn_ratio_mean(rails_dur,:)=y2;
        turn_ratio_std(rails_dur,:)=sem;
        turn_ratio_num_worms(rails_dur,:)=n_tracks;
    elseif mn==3
        reverse_ratio_mean(rails_dur,:)=y2;
        reverse_ratio_std(rails_dur,:)=sem;
        reverse_ratio_num_worms(rails_dur,:)=n_tracks;
    end
    
end

end

end  %%% this will prevent plotting any graphs and the corresponding calculations
end

%%
if plot_comparison_reverse_beh==1
    
my_colors = lines(size(rails_durations,2));
figure;
for rails_dur=1:size(rails_durations,2)
    errorbar(stimulus_intensity_array, reverse_ratio_mean(rails_dur,:),reverse_ratio_std(rails_dur,:), 'o-', 'color', my_colors(rails_dur,:),'LineWidth',2,'Markersize',10)
    hold on
    ax = gca;
    ax.XTick = stimulus_intensities;
    % ax.YTick = [0 0.3 0.6];
    ax.FontSize = 13;

    xlabel('Stimulus Intensity (\muW/mm^2)') % x-axis label
    ylabel([behavior_names{3}, ' Ratio'])
    title(['Beh. = ', behavior_names{3},...
            ',' ' Stim width = ', num2str(current_rails_dur/parameters.SampleRate),' s']);
%     legend({[num2str(rails_durations(1)), 'sec'],[num2str(rails_durations(2)), 'sec'],[num2str(rails_durations(3)), 'sec']},'Location','northwest')
%     ylim([0,0.7]);
    grid off;
    xlim([0 max(stimulus_intensities)+2])
    box off
end
hold off;

figure;
for rails_dur=1:size(rails_durations,2)

    errorbar(stimulus_intensity_array, turn_ratio_mean(rails_dur,:),turn_ratio_std(rails_dur,:), 'o-', 'color', my_colors(rails_dur,:),'LineWidth',2,'Markersize',10)
    hold on
    ax = gca;
    ax.XTick = stimulus_intensities;
    % ax.YTick = [0 0.3 0.6];
    ax.FontSize = 13;

    xlabel('Stimulus Intensity (\muW/mm^2)') % x-axis label
    ylabel([behavior_names{2}, ' Ratio'])
    title(['Beh. = ', behavior_names{2},...
            ',' ' Stim width = ', num2str(current_rails_dur/parameters.SampleRate),' s']);
%     legend({[num2str(rails_durations(1)), 'sec'],[num2str(rails_durations(2)), 'sec'],[num2str(rails_durations(3)), 'sec']},'Location','northwest')
    grid off;
%     ylim([0.7,1]);
    xlim([0 max(stimulus_intensities)+2])
    box off

end
hold off;
end

if save_mat_file==1
    if analysis_rails==1
    writematrix(count_number_of_reversal_folder_rails,excel_filename);
    elseif analysis_stim_while_turning==1
    writematrix(count_number_of_reversal_folder_turns,excel_filename);
    end  
end

%%
%%%%% coming up with a way, so that I can use my behavior annotation
for ijk=1:size(all_behavior_annotations_for_frame,2)
    for mn=1:size(all_behavior_annotations_for_frame{1,ijk},2)
        for pq=1:size(all_behavior_annotations_for_frame{1,ijk}{1,mn},2)
            sk_behavior_annotations_for_frame{1,ijk}{1,mn}(1,pq)=1; 
            if all_velocity_annotations_for_frame{1,ijk}{1,mn}(1,pq)<0
                sk_behavior_annotations_for_frame{1,ijk}{1,mn}(1,pq)=3;
            end 
            
% % % % %             if all_velocity_annotations_for_frame{1,ijk}{1,mn}(1,pq)<0.02 && all_velocity_annotations_for_frame{1,ijk}{1,mn}(1,pq)>-0.02
% % % % %                 sk_behavior_annotations_for_frame{1,ijk}{1,mn}(1,pq)=4;
% % % % %             end

            
            if all_ellipse_ratio_annotations_for_frame{1,ijk}{1,mn}(1,pq)<ellipse_ratio_threshold
                sk_behavior_annotations_for_frame{1,ijk}{1,mn}(1,pq)=2;
            end 
        end
    end
end

original_behavior_annotations_for_frame=all_behavior_annotations_for_frame;
all_behavior_annotations_for_frame=sk_behavior_annotations_for_frame;

%%
if plot_behavior_heatmaps==1 && plotting_graphs==1 
    
    array_to_withdraw_random_samples_from=[];
    for m=[1 size(all_behavior_annotations_for_frame,2)]
        array_to_withdraw_random_samples_from(m,1)=size(all_behavior_annotations_for_frame{1,m}{1,1},2);
    end
    array_to_withdraw_random_samples_from_final=min(nonzeros(array_to_withdraw_random_samples_from));
           
    color_matrix=behavior_colors;
    color_matrix(1,:)=[0.9 0.9 0.9];
    f = figure;
    f.Position = [100 300 1140 520];
 
    for ij=1:size(stimulus_intensities,2)

        h3=subplot(1,size(stimulus_intensities,2),ij);
% % % %         %%% skipping plotting for 40uW
% % % %         if ij==2
% % % %             ij=3;
% % % %         end
        
        hAxes = gca;
        clims = [0 3.0001];
        dummy_data_for_heatmap=all_behavior_annotations_for_frame{1,ij}(:,1:(2*current_rails_dur)+1);
        behavior_data_for_heatmap=cell2mat(reshape(dummy_data_for_heatmap,[(2*current_rails_dur+1),1]));
        behavior_data_for_heatmap=behavior_data_for_heatmap';
        behavior_data_for_heatmap_cell{1,ij}=behavior_data_for_heatmap;
        behavior_data_for_heatmap_cell{2,ij}=stimulus_intensities(1,ij);
        
        if analysis_rails==1        
            rng(random_seed_open_loop_heatmap,'twister')
        else
            rng(random_seed_closed_loop_heatmap,'twister')
        end
        random_samples=randsample(array_to_withdraw_random_samples_from_final,number_of_rastors);
        
        data_to_be_plotted=behavior_data_for_heatmap(random_samples,:);
        [sorted_value,sorted_idx] = sort_heatmap_data(data_to_be_plotted,round(size(data_to_be_plotted,2)/2),(round(size(data_to_be_plotted,2)/2)+current_rails_dur),behavior_sort_threshold,'ascend',behavior_comparision_mode);
        imagesc(hAxes,data_to_be_plotted(sorted_idx,:),clims);
%         imagesc(hAxes,data_to_be_plotted,clims);
        hold on;
        xline(round(size(data_to_be_plotted,2)/2),'--k','linewidth',1.0);
        hold on;
        xline(round(size(data_to_be_plotted,2)/2)+current_rails_dur,'--k','linewidth',1.0);
        colormap( hAxes , color_matrix)
        yticks([10:10:30])
        
        xticks([1:(floor(size(data_to_be_plotted,2)/2)):size(data_to_be_plotted,2)])
        xt = get(gca, 'XTick');                                             % Original 'XTick' Values
        xtlbl = linspace(-rails_dur_itr, rails_dur_itr, numel(xt));                     % New 'XTickLabel' Vector
        set(gca, 'XTick',xt, 'XTickLabel',xtlbl)   % Label Ticks
        
        xt = get(gca, 'XTick');                                             % Original 'XTick' Values
        xtlbl = linspace(-5, 5, numel(xt));                     % New 'XTickLabel' Vector
        set(gca, 'XTick',xt, 'XTickLabel',xtlbl)   % Label Ticks
        title([num2str(stimulus_intensities(ij)) '\muW/mm^2'])
        ylabel('Worm #')
        xlabel('Time (s)')
        box off
        ax = gca; ax.FontSize = 20;
        
        % Create rectangle
        annotation(f,'rectangle', [0.297 0.94 0.101 0.02], 'Color',[0 0 0], 'FaceColor',[1 1 1],'LineWidth',0.5);
        annotation(f,'rectangle', [0.737 0.94 0.101 0.02], 'Color',[0 0 0], 'FaceColor',[1 0 0],'LineWidth',0.5);

        if ij ~= 1
           set(gca,'ytick',[])
           set(gca,'yticklabel',[])
           set(gca,'ylabel',[])
           set(gca,'xlabel',[])
        end
        if ij==3
            cb=colorbar('position',[.93 .18 .02 .745],'Ticks',[0 1 2 3],'FontSize',16);
%             originalSize2 = get(gca, 'Position');
%             set(h3, 'Position', originalSize2); 
        end
    end
end

%% plotting velocity plots

if plot_velocity_traces==1 && plotting_graphs==1
    f = figure;
    f.Position = [100 300 1140 520];
    time_axis=[(-1*rails_dur_itr):1/fps:rails_dur_itr];
    n_bins = 60;
    edges = linspace(-0.3,0.3,n_bins);
    velocity_density = zeros(numel(edges)-1, numel(time_axis));
    for ij=1:size(stimulus_intensities,2)
        h3=subplot(1,size(stimulus_intensities,2),ij);
% %         if ij==2
% %             ij=3;
% %         end
       
       dummy_velocity_data_for_heatmap=all_velocity_annotations_for_frame{1,ij}(:,1:(2*current_rails_dur)+1);
       behavior_data_for_heatmap=cell2mat(reshape(dummy_velocity_data_for_heatmap,[(2*current_rails_dur+1),1]));
       behavior_data_for_heatmap=behavior_data_for_heatmap';
       mean_velocity=mean(behavior_data_for_heatmap,1);  
       
       for time_index = 1:length(time_axis)
           velocity_density(:,time_index) = histcounts(behavior_data_for_heatmap(:,time_index), edges,'Normalization','probability');
       end
        
       hold on
       imagesc(time_axis, edges, velocity_density);
       hold on;
       plot(time_axis,mean_velocity,'-k','LineWidth',2);
       hold on;
       xline(0,'--y','linewidth',1);
       hold on;
       xline(rails_dur_itr,'--y','linewidth',1);
       hold on;
       title([num2str(stimulus_intensities(ij)) '\muW/mm^2'])
       ylim([-0.3 0.3])
       xlim([-rails_dur_itr rails_dur_itr])
       yticks([-0.3:0.3:0.3])
       xticks([-rails_dur_itr:rails_dur_itr:rails_dur_itr])
       ylabel('Velocity (mm/sec)')
       xlabel('Time (s)')
% %        cb=colorbar;
       box off
       ax = gca; ax.FontSize = 20;
        annotation(f,'rectangle', [0.297 0.94 0.101 0.02], 'Color',[0 0 0], 'FaceColor',[1 1 1],'LineWidth',0.5);
        annotation(f,'rectangle', [0.737 0.94 0.101 0.02], 'Color',[0 0 0], 'FaceColor',[1 0 0],'LineWidth',0.5);
       if ij ~= 1
            set(gca,'ytick',[])
            set(gca,'yticklabel',[])
            set(gca,'ylabel',[])
            set(gca,'xlabel',[])
%             set(gca,'xtick',[])
%             set(gca,'xticklabel',[])
       end
       if ij==3
            cb=colorbar('position',[.93 .18 .02 .745],'Ticks',[0 0.1],'FontSize',16);
% %             originalSize3 = get(gca, 'Position');
% %             set(h3, 'Position', originalSize3); 
       end
    end
end

%%
%%%% plotting ellipse ratio traces
if plot_ellipse_ratio_traces==1 && plotting_graphs==1
    f = figure;
    f.Position = [100 300 1140 520];
    time_axis=[(-1*rails_dur_itr):1/fps:rails_dur_itr];
    n_bins = 60;
    er_edges = linspace(1,7,n_bins);
    er_density = zeros(numel(er_edges)-1, numel(time_axis));
    for ij=1:size(stimulus_intensities,2)
        h3=subplot(1,size(stimulus_intensities,2),ij);
% % %         if ij==2
% % %             ij=3;
% % %         end
       
       dummy_ellipse_ratio_data_for_heatmap=all_ellipse_ratio_annotations_for_frame{1,ij}(:,1:(2*current_rails_dur)+1);
       ellipse_ratio_data_for_heatmap=cell2mat(reshape(dummy_ellipse_ratio_data_for_heatmap,[(2*current_rails_dur+1),1]));
       ellipse_ratio_data_for_heatmap=ellipse_ratio_data_for_heatmap';
       mean_ellipse_ratio=mean(ellipse_ratio_data_for_heatmap,1);  

       for time_index = 1:length(time_axis)
           er_density(:,time_index) = histcounts(ellipse_ratio_data_for_heatmap(:,time_index), er_edges,'Normalization','probability');
       end
        
       hold on
       imagesc(time_axis, er_edges, er_density);
       hold on;
       plot(time_axis,mean_ellipse_ratio,'-k','LineWidth',2);
       hold on;
       xline(0,'--y','linewidth',1);
       hold on;
       xline(rails_dur_itr,'--y','linewidth',1);
       hold on;
% %        title([num2str(stimulus_intensities(ij)) '\muW/mm^2'])
       ylim([1 7])
       xlim([-rails_dur_itr rails_dur_itr])
       yticks([2:2:6])
       xticks([-rails_dur_itr:rails_dur_itr:rails_dur_itr])
       ylabel('Ellipse Ratio')
       xlabel('Time (s)')
% %        cb=colorbar;
       box off
       ax = gca; ax.FontSize = 20;
       annotation(f,'rectangle', [0.297 0.94 0.101 0.02], 'Color',[0 0 0], 'FaceColor',[1 1 1],'LineWidth',0.5);
       annotation(f,'rectangle', [0.737 0.94 0.101 0.02], 'Color',[0 0 0], 'FaceColor',[1 0 0],'LineWidth',0.5);
       if ij ~= 1
            set(gca,'ytick',[])
            set(gca,'yticklabel',[])
            set(gca,'ylabel',[])
            set(gca,'xlabel',[])
%             set(gca,'xtick',[])
%             set(gca,'xticklabel',[])
       end
       if ij==3
            cb=colorbar('position',[.93 .18 .02 .745],'Ticks',[0 0.1],'FontSize',16);
%             originalSize3 = get(gca, 'Position');
%             set(h3, 'Position', originalSize3); 
       end
    end
end

%%
if plot_fwd_to_turn_heatmaps==1 && plotting_graphs==1 && analysis_stim_while_turning==1

color_matrix=behavior_colors;
color_matrix(1,:)=[1 1 1];
f = figure;
f.Position = [100 300 1140 520];
for ij=1:size(stimulus_intensities,2)

    h3=subplot(1,size(stimulus_intensities,2),ij);

% % % %     %%% skipping plotting for 40uW
% % % %     if ij==2
% % % %         ij=3;
% % % %     end

    hAxes = gca;
    clims = [0 3.001];
    dummy_data_for_heatmap=all_behavior_annotations_for_frame{1,ij}(:,1:(2*current_rails_dur)+1);
    behavior_data_for_heatmap=cell2mat(reshape(dummy_data_for_heatmap,[(2*current_rails_dur)+1,1]));
    behavior_data_for_heatmap=behavior_data_for_heatmap';
    
    fwd_to_turn_data=[];
    for i=1:size(behavior_data_for_heatmap,1)
        if any(behavior_data_for_heatmap(i,61:120)~=1)
            continue
        end
        fwd_to_turn_data=[fwd_to_turn_data; behavior_data_for_heatmap(i,:)];
    end

    rng(random_seed_forward_to_turn_heatmap,'twister')
    random_samples=randsample(size(fwd_to_turn_data,1),number_of_rastors);
    data_to_be_plotted=fwd_to_turn_data(random_samples,:);
    [sorted_value,sorted_idx] = sort_heatmap_data(data_to_be_plotted,round(size(data_to_be_plotted,2)/2),(round(size(data_to_be_plotted,2)/2)+current_rails_dur),behavior_sort_threshold,'ascend',behavior_comparision_mode);
    imagesc(hAxes,data_to_be_plotted(sorted_idx,:),clims);

    hold on;
    xline(151,'--k','linewidth',1);
    hold on;
    xline(241,'--k','linewidth',1);
    colormap( hAxes , color_matrix)
    xticks([1:floor(total_window_for_visualization/2):total_window_for_visualization])
    xt = get(gca, 'XTick');                                             % Original 'XTick' Values
    xtlbl = linspace(-total_time_for_visualization, total_time_for_visualization, numel(xt));                     % New 'XTickLabel' Vector
    set(gca, 'XTick',xt, 'XTickLabel',xtlbl)   % Label Ticks
    xt = get(gca, 'XTick');                                             % Original 'XTick' Values
    xtlbl = linspace(-5, 5, numel(xt));                     % New 'XTickLabel' Vector
    set(gca, 'XTick',xt, 'XTickLabel',xtlbl)   % Label Ticks
    yticks([10:10:30])
%     title([num2str(stimulus_intensities(ij)) '\muW/mm^2'])
    ylabel('Worm #')
    xlabel('Time (s)')
    box off
    ax = gca; ax.FontSize = 20;
    annotation(f,'rectangle',[0.297 0.94 0.101 0.02],'Color',[0 0 0],'FaceColor',[1 1 1],'LineWidth',0.5);
    annotation(f,'rectangle',[0.737 0.94 0.101 0.02],'Color',[0 0 0],'FaceColor',[1 0 0],'LineWidth',0.5);
    if ij ~= 1
       set(gca,'ytick',[])
       set(gca,'yticklabel',[])
       set(gca,'ylabel',[])
       set(gca,'xlabel',[])
    end
    if ij==3
        cb=colorbar('position',[.93 .18 .02 .745],'Ticks',[0 1 2 3]);
%         originalSize2 = get(gca, 'Position');
%         set(h3, 'Position', originalSize2); 
    end
end
end

%%

if plot_reverse_to_turn_heatmaps==1 && plotting_graphs==1 && analysis_stim_while_turning==1

color_matrix=behavior_colors;
color_matrix(1,:)=[0.9 0.9 0.9];
f = figure;
f.Position = [100 300 1140 520]; %%%%[100 100 1140 520];
for ij=1:size(stimulus_intensities,2)

    h3=subplot(1,size(stimulus_intensities,2),ij);

% % %     %%% skipping plotting for 40uW
% % %     if ij==2
% % %         ij=3;
% % %     end

    hAxes = gca;
    clims = [0 3.001];
    dummy_data_for_heatmap=all_behavior_annotations_for_frame{1,ij};
    behavior_data_for_heatmap=cell2mat(reshape(dummy_data_for_heatmap,[size(dummy_data_for_heatmap,2),1]));
    behavior_data_for_heatmap=behavior_data_for_heatmap';
    
    reverse_to_turn_data=[];
    for i=1:size(behavior_data_for_heatmap,1)
        if sum(behavior_data_for_heatmap(i,31:120)==3)<=60
            continue
        end
        reverse_to_turn_data=[reverse_to_turn_data; behavior_data_for_heatmap(i,:)];
    end

    rng(random_seed_reverse_to_turn_heatmap,'twister')
    
    random_samples=randsample(size(reverse_to_turn_data,1),number_of_rastors);
    data_to_be_plotted=reverse_to_turn_data(random_samples,:);
    [sorted_value,sorted_idx] = sort_heatmap_data(data_to_be_plotted,round(size(data_to_be_plotted,2)/2),(round(size(data_to_be_plotted,2)/2)+current_rails_dur),behavior_sort_threshold,'ascend',behavior_comparision_mode);
    imagesc(hAxes,data_to_be_plotted(sorted_idx,:),clims);

    hold on;
    xline(151,'--k','linewidth',1);
    hold on;
    xline(241,'--k','linewidth',1);
    colormap( hAxes , color_matrix)
    xticks([1:floor(total_window_for_visualization/2):total_window_for_visualization])
    xt = get(gca, 'XTick');                                             % Original 'XTick' Values
    xtlbl = linspace(-total_time_for_visualization, total_time_for_visualization, numel(xt));                     % New 'XTickLabel' Vector
    set(gca, 'XTick',xt, 'XTickLabel',xtlbl)   % Label Ticks

    xt = get(gca, 'XTick');                                             % Original 'XTick' Values
    xtlbl = linspace(-5, 5, numel(xt));                     % New 'XTickLabel' Vector
    set(gca, 'XTick',xt, 'XTickLabel',xtlbl)   % Label Ticks
    yticks([10:10:30])
%     title([num2str(stimulus_intensities(ij)) '\muW/mm^2'])
    ylabel('Worm #')
    xlabel('Time (s)')
    box off
   ax = gca; ax.FontSize = 18;
   annotation(f,'rectangle',[0.297 0.94 0.101 0.02],'Color',[0 0 0],'FaceColor',[1 1 1],'LineWidth',0.5);
   annotation(f,'rectangle',[0.737 0.94 0.101 0.02],'Color',[0 0 0],'FaceColor',[1 0 0],'LineWidth',0.5);
    if ij ~= 1
       set(gca,'ytick',[])
       set(gca,'yticklabel',[])
       set(gca,'ylabel',[])
       set(gca,'xlabel',[])
    end
    if ij==3
        cb=colorbar('position',[.93 .18 .02 .745],'Ticks',[0 1 2 3]);
% %         originalSize2 = get(gca, 'Position');
% %         set(h3, 'Position', originalSize2); 
    end
end
end

%% to count the probability of reversals
if count_prob_reversal==1
   for ij=1:size(stimulus_intensities,2)
       
% % %     %%% skipping plotting for 40uW
% % %     if ij==2
% % %         ij=3;
% % %     end

    dummy_data_to_count_prob_reversal=all_velocity_annotations_for_frame{1,ij}(:,1:(2*current_rails_dur)+1);
    velocity_data_to_count_prob_reversal=cell2mat(reshape(dummy_data_to_count_prob_reversal,[(2*current_rails_dur)+1,1]));
    velocity_data_to_count_prob_reversal=velocity_data_to_count_prob_reversal';

    for uv=1:size(velocity_data_to_count_prob_reversal,1)
        if any(velocity_data_to_count_prob_reversal(uv,round(size(velocity_data_to_count_prob_reversal,2)/2):(round(size(velocity_data_to_count_prob_reversal,2)/2)+current_rails_dur))<=negative_vel_threshold)
            count_number_of_reversals_array{1,ij}(uv,1)=1;
        else
            count_number_of_reversals_array{1,ij}(uv,1)=0;
        end
    end
    
    [mean_data(:,ij),err_high(:,ij),err_low(:,ij)]=bootstrap_mean_and_ci(10000,0.05,count_number_of_reversals_array{1,ij});
    
   end 
   
   number_of_reversals_and_stim_intensity=[count_number_of_reversals_array; num2cell(stimulus_intensities)];
   
end

end
