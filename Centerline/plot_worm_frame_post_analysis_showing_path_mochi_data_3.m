function [] = plot_worm_frame_post_analysis_showing_path_mochi_data_3(I, center_lines, centerline_color, centerline_stimulus, debugimage,in_track_index,current_ellipse_ratio,smoothed_path_x,smoothed_path_y,stimulus_delivered,path_history,demo_mode,all_velocity,current_rails_dur,initial_frame,stim_color)

% % % % % % [] = plot_worm_frame(I, center_lines, centerline_color, UncertainTips, centerline_stimulus, debugimage,in_track_index,current_ellipse_ratio)

%     IWFig = findobj('Tag', ['IWFig', num2str(plotting_index)]);
%     if isempty(IWFig)
%         IWFig = figure('Tag', ['IWFig', num2str(plotting_index)]);
%     else
%         figure(IWFig);
%     end
    %used for debugging
%     hold off;

%     clf

    f(1)=subplot(1,3,[1 2]);
    
    pos1 = [-0.1 0.1 0.65 0.65];
    subplot('Position',pos1)
    
    if nargin > 5
        
        %debugimage inputted
        I=imadjust(I);
        I = I ./ 1.1;
        imshow(cat(3, I, I, I) + debugimage, [], 'InitialMagnification', 500, 'Border','tight');
        
%         imshow(I, [], 'InitialMagnification', 500, 'Border','tight'); hold on;
% %         if ~isempty(find(debugimage(:,:,1)>10))
% %             hold on; 
% % %             text(60, 0, num2str('.'), 'Color', 'r','FontSize',140);
% %             LEDPower=90;
% %             [frame_h, frame_w] = size(I);
% %             plot_x = ceil(frame_w - (frame_w/10));
% %             plot_y = ceil(frame_h/10);
% %             plot(plot_x-5, plot_y+4, 'o', 'MarkerSize', 30, 'MarkerEdgeColor','none', 'MarkerFaceColor',[1 0 0]./1)
% %             text(73,4, [num2str(round(LEDPower)), ' uW/mm^{2}'], 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', 'color', 'y', 'fontsize', 12);
% %         end
    else
        imshow(I, [], 'InitialMagnification', 500, 'Border','tight');
    end
    
    %%% path
    smoothed_path=[smoothed_path_x smoothed_path_y];
    
    if isempty(centerline_color)
        centerline_color = [1, 0, 0];
    end
    
%     center_lines=flip(center_lines);
    centerline_correction=fliplr(center_lines(10,:));
    
    %%%% plot the previous few seconds of path or else plot from the beginning
    if in_track_index>path_history
        path_t = smoothed_path(in_track_index-path_history:in_track_index,:) - smoothed_path(in_track_index,:)+centerline_correction;
    else
        path_t = smoothed_path(1:in_track_index,:) - smoothed_path(in_track_index,:)+centerline_correction;
    end

    hold on
% % %     %head
    plot(center_lines(1,2), center_lines(1,1), '.', 'Color', centerline_color, 'markersize',40)
    plot(center_lines(:,2), center_lines(:,1), '-', 'Color', centerline_color, 'LineWidth',4)
% %     if ~isempty(centerline_stimulus)
% %         for centerline_point_index = 1:size(center_lines,1)
% %             if isnan(centerline_stimulus(centerline_point_index))
% %                 dot_color = [0 1 0];
% %             elseif centerline_stimulus(centerline_point_index) > 0
% %                 dot_color = [centerline_stimulus(centerline_point_index)/255, 0, 0];
% %             elseif centerline_stimulus(centerline_point_index) < 0
% %                 dot_color = [0,0,-centerline_stimulus(centerline_point_index)/255];
% %             else
% %                 dot_color = [0 0 0];
% %             end
% %             plot(center_lines(centerline_point_index,2), center_lines(centerline_point_index,1), '.', 'Color', dot_color, 'markersize',20)
% %         end
% %     end
    hold on;
    plot(path_t(:,1),path_t(:,2),'.-g','LineWidth',4);
    hold on;
    
    plot(path_t(1:30:end,1),path_t(1:30:end,2),'.', 'MarkerSize',18,'color',[1 1 0],'LineWidth',2);
    hold on;
        

% % % % %     %uncertain tips
% % % % %     if ~isempty(UncertainTips) && ~isempty(UncertainTips.Tips)
% % % % %         plot(UncertainTips.Tips(:,2), UncertainTips.Tips(:,1), 'oy')
% % % % %     end
%     %tail
%     plot(center_lines(end,2), center_lines(end,1), 'ob')
% 
%     
% %     %direction
% % %     quiver(size(I,2)/2, size(I,1)/2, sind(direction)*speed*100, -cosd(direction)*speed*100, 'AutoScale','off', 'Linewidth', 1.5);
% % %     %title (['Eccentricity = ', num2str(eccentricity)]);    
% % %     title (sprintf('Frame=%d, E.R.=%.1f',in_track_index,current_ellipse_ratio));
    
    if demo_mode~=1
    %score
    text(10, 10, num2str(in_track_index), 'Color', 'y','FontSize',14);
    %eccentricity
    text(10, 60, num2str(current_ellipse_ratio,'%4.2f'), 'Color', 'y','FontSize',14);
    
    %%% stim intensity
% %     text(60,10, [num2str(round(abs(stimulus_delivered))), ' uW/mm^{2}'], 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', 'color', 'y', 'fontsize', 14);

    hold off
    end
    hold off;
    
    f(2)=subplot(1,3,3);    
    pos2 = [0.55 0.2 0.35 0.55];
    subplot('Position',pos2);

    plot(all_velocity,'-r','Linewidth',2)
% % %     plot([-round(size(all_velocity,1)./2)-1:1:round(size(all_velocity,1)./2)],all_velocity,'-r','Linewidth',2)
    xlim([initial_frame-120 initial_frame+120])
    ylim([-0.4 0.4])
    hold on;
    xline(in_track_index,'color','b','Linewidth',2);
    hold on;
    yline(0,'color','k','Linewidth',1);
    xlabel('Frames','FontSize',14)
    ylabel('Velocity (mm/s)','FontSize',14)
    hold on;
    if in_track_index>=initial_frame && in_track_index-initial_frame<=current_rails_dur
        light_pulse = area([initial_frame initial_frame+(in_track_index-initial_frame)], [7 7],'BaseValue',-0.4,'ShowBaseLine','off','FaceColor',stim_color,'LineStyle','none');
        alpha(0.2);
    hold on
    end
    
    if in_track_index>initial_frame+current_rails_dur
        light_pulse = area([initial_frame (initial_frame+current_rails_dur)], [7 7],'BaseValue',-0.4,'ShowBaseLine','off','FaceColor',stim_color,'LineStyle','none');
        alpha(0.2);
        hold on
    end
    
    ax = gca;
    ax.FontSize = 14;
    hold on;
    set(gcf,'color','w');
%     hold off;
    
%     drawnow
end