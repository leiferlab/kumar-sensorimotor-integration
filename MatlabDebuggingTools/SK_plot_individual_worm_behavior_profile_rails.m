clear
clc
close all

main_folder=('/projects/LEIFER/Sandeep/APIData/20200831_RunRailsTriggeredByTurning_Sandeep_AML17_AVA_10ulRet/Data20200831_162748'); %% main folder path
filename='workspace_data_20200831_162748_80uW_5sec_turning.mat';
load(fullfile(main_folder,filename)); %% save the workspace data inside the main folder

%%%%%%%%% user defined inputs
plot_figures=11;
to_debug_wrong_trials=1;


if contains(filename,'full')
    testtype='rails';
elseif contains(filename,'turn')
    testtype='turns';
end

if testtype=='rails'
rails_indices=fullrails_peaks_array_excel_data;
% rails_indices=cell2mat(fullrails_peaks_array_excel_data);
elseif testtype=='turns'
rails_indices=cell2mat(stim_while_turning_peaks_array_excel_data);
end

for mn=1:size(rails_indices,1)
    uniqueNum = unique(rails_indices(mn,:));
    [n,bin] = histc(rails_indices(mn,:),uniqueNum);
    if any(n(2:end)>1)
        dummy_rails_indices(mn,:)=[uniqueNum zeros(size(rails_indices,2)-size(uniqueNum,2))];
    else
        dummy_rails_indices(mn,:)=rails_indices(mn,:);
    end
end

rails_indices=dummy_rails_indices;
%%
for i=1:size(rails_indices,1)
% for i=[5]
%     i
    for j=2:size(rails_indices,2)
%     for j=2
    if rails_indices(i,j)==0
        continue
    end
    
    if size(current_tracks(rails_indices(i,1)).BehavioralAnnotation,2) < rails_indices(i,j)+500
        sprintf('trial data size is insufficient')
        continue 
    end
    
    rails_indices_nonzero(i,:)=rails_indices(i,:);
    
    original_data=current_tracks(rails_indices(i,1)).BehavioralAnnotation(1,(rails_indices(i,j)-500):(rails_indices(i,j)+500));
    original_data_matrix{i,j-1}=original_data;
        
    original_velocity=current_tracks(rails_indices(i,1)).Velocity((rails_indices(i,j)-500):(rails_indices(i,j)+500),1);
    original_velocity=smooth(original_velocity,5);
    original_velocity=original_velocity';
    original_velocity_matrix{i,j-1}=original_velocity;
   
    original_direction=current_tracks(rails_indices(i,1)).Direction(1,(rails_indices(i,j)-500):(rails_indices(i,j)+500));
    
    data_reverse=find(original_data>2.5); %%%% when worm is reversing
        
    data_reverse_final=zeros(1,1000);
    data_reverse_final(data_reverse)=3;
        
    [~, reverse_locs,reverse_widths,~] = findpeaks(data_reverse_final);
    
    %%%%%% if the worm is reversing continuously for more than six seconds
    %%%%%% then reassign those elements as movinf forward
    if any(reverse_widths>270)
        index=find(reverse_widths>270);
        reverse_locs_itr=reverse_locs(index);
        reverse_widths_itr=reverse_widths(index);
        original_data(reverse_locs_itr:reverse_locs_itr+reverse_widths_itr)=1;
    end
    
    %%%%if pipeline does'nt assign behavior to one of the three states then we will take velocity into consideration to assign forward or reverse
    data_unassigned=find(original_data==0); %%%% when beh is unassigned

    corrected_data=original_data;
    for xx=[data_unassigned]
        if original_velocity(1,xx)<=-0.05
            corrected_data(1,xx)=3;
        elseif original_velocity(1,xx)>=0.1
            corrected_data(1,xx)=1;
        elseif original_velocity(1,xx)>-0.05 && original_velocity(1,xx)<0.1
            corrected_data(1,xx)=2;
        end
    end
    
    corrected_data_matrix{i,j-1}=corrected_data;
        
    if plot_figures==1
        figure;
        x_data=(rails_indices(i,j)-500):1:(rails_indices(i,j)+500);
        plot(x_data,original_data,'--g')
        hold on;
        plot(x_data,corrected_data,'-r')
        hold on;
        plot(x_data,original_velocity,'-b')
        hold on;
%         plot(x_data,original_direction/60,'-g')
%         hold on;
        title(['worm ID: ' num2str(rails_indices(i,1)) ',' ' track id: ' num2str(rails_indices(i,j))])
        xlim([rails_indices(i,j)-500 rails_indices(i,j)+500])
        xlabel('Time (sec)'); ylabel('Behavior annotation');
        xline(rails_indices(i,j),'--k');
        hold on;
        xline(rails_indices(i,j)+150,'--k');
        hold on;
        yline(0,'--k');
        hold off;
        ax = gca;
        ax.FontSize = 14;
        box off;
        axis tight;
    end
   
    behavior_transitions_table{i,1}=rails_indices(i,1);
    behavior_transitions_table{i,2*j}=rails_indices(i,j);
    
    if max(original_velocity)<=0.16 && min(original_velocity)>=-0.15  %%%% ignore trials when worms hardly moved
        behavior_transitions_table{i,(2*j+1)}='i,';
        ignore_type_matrix{i,j}='type 1';
       continue
    end
    
    if max(original_velocity(:,1:400))<=0.06 && min(original_velocity(:,1:400))>=-0.06  %%%% ignore trials when worms hardly moved
        behavior_transitions_table{i,(2*j+1)}='i,';
        ignore_type_matrix{i,j}='type 2';
       continue
    end
        
    if testtype=='rails'  %%%% ignore these trials only when the stim happens at reversal or turns
    if all(corrected_data(480:510)==2) || all(corrected_data(480:510)==3)
        behavior_transitions_table{i,(2*j+1)}='i,';
        ignore_type_matrix{i,j}='type 3';
       continue
    end
    end
    
%     if testtype=='turns'  %%%% ignore these trials only when the stim happens at reversal 
%     if all(corrected_data(501:530)==3)
%         behavior_transitions_table{i,(2*j+1)}='i,';
%         ignore_type_matrix{i,j}='type 4';
%        continue
%     end
%     end
%     
    if testtype=='rails'
        if any(corrected_data(511:560)==3) %%%% maybe make it 511:600
           behavior_transitions_table{i,(2*j+1)}='t,';
%         elseif any(corrected_data(501:650)==2)
%            behavior_transitions_table{i,(2*j+1)}='t,';
        else
            behavior_transitions_table{i,(2*j+1)}='f,';
        end
        
    elseif testtype=='turns'
        if any(corrected_data(525:600)==3)
            behavior_transitions_table{i,(2*j+1)}='t,';
        else
        behavior_transitions_table{i,(2*j+1)}='f,';
        end
    end
        
    end
end

behavior_transitions_table(all(cellfun('isempty',behavior_transitions_table),2),:) = [];
behavior_transitions_table(:,all(cellfun('isempty',behavior_transitions_table),1),:) = [];

corrected_data_matrix(all(cellfun('isempty',corrected_data_matrix),2),:) = [];
corrected_data_matrix(:,all(cellfun('isempty',corrected_data_matrix),1),:) = [];

original_data_matrix(all(cellfun('isempty',original_data_matrix),2),:) = [];
original_data_matrix(:,all(cellfun('isempty',original_data_matrix),1),:) = [];

original_velocity_matrix(all(cellfun('isempty',original_velocity_matrix),2),:) = [];
original_velocity_matrix(:,all(cellfun('isempty',original_velocity_matrix),1),:) = [];

rails_indices_final = rails_indices_nonzero(any(rails_indices_nonzero,2),:);
%%
%%%%%%%% to plot figures for the incorrect precticted trials
if to_debug_wrong_trials==1
    
    [worm_ids_dummy,rails_150uW_5sec,~]=xlsread('AML17_data_20200831_162748_80uW_5sec_reanalyzed_debug.xlsx');
    true_data=rails_150uW_5sec(2:end,3:2:end);
    computed_data=behavior_transitions_table(:,3:2:end);
    to_compare_data=[true_data computed_data];
%%
 
    for i=1:size(to_compare_data,1)
%     for i=[1 2 3 4 7 11]
        for j=1:(size(to_compare_data,2)/2)
            if size(to_compare_data{i,j},1)==0
                continue
            end
        
        if to_compare_data{i,j}==to_compare_data{i,(size(to_compare_data,2)/2)+j}
            to_compare_table_results{i,j}=1;
        else
            to_compare_table_results{i,j}=0;
        end
        end
    end

    empties_table_results = cellfun('isempty',to_compare_table_results);
    to_compare_table_results(empties_table_results) = {NaN};
    to_compare_table_results=cell2mat(to_compare_table_results);
    correct_annotation=sum(to_compare_table_results(:) == 1);
    wrong_annotation=sum(to_compare_table_results(:) == 0);

    index_wrong_annotation=find(to_compare_table_results==0)';
    num_worms_analyzed=sum(sum(~isnan(to_compare_table_results)));
    percent_correct_analysis=(num_worms_analyzed-size(index_wrong_annotation,2))/num_worms_analyzed*100
    
    compiled_comparision_results=[num2cell(rails_indices_final(:,1)) to_compare_data];
    rails_indices_debug=rails_indices_final(:,2:end);

for xx=1:size(to_compare_table_results,1)
    for yy=1:size(to_compare_table_results,2)
        if to_compare_table_results(xx,yy)==0
        corrected_data_debug=corrected_data_matrix{xx,yy};
        original_velocity_debug=original_velocity_matrix{xx,yy};
        original_data_debug=original_data_matrix{xx,yy};
%         xx
%         yy
        figure;
        x_data=(rails_indices_debug(xx,yy)-500):1:(rails_indices_debug(xx,yy)+500);
        plot(x_data,corrected_data_debug,'-r')
        hold on;
        plot(x_data,original_data_debug,'--g')
        hold on;
        plot(x_data,original_velocity_debug,'-b')
        hold on;
        title(['worm ID: ' num2str(rails_indices_final(xx,1)) ',' ' track id: ' num2str(rails_indices_debug(xx,yy))])
        xlim([rails_indices_debug(xx,yy)-500 rails_indices_debug(xx,yy)+500])
        xlabel('Time (sec)'); ylabel('Behavior annotation');
        xline(rails_indices_debug(xx,yy),'--k');
        hold on;
        xline(rails_indices_debug(xx,yy)+150,'--k');
        hold on;
        yline(0,'--k');
        hold off;
        ax = gca;
        ax.FontSize = 14;
        box off;
        axis tight;
        continue
        end
    end
end
end