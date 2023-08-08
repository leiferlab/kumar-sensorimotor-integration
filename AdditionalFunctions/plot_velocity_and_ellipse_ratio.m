function plot_velocity_and_ellipse_ratio(stim_while_turning_peaks,current_rails_dur,dummy_array_with_velocity_info,dummy_array_with_ellipse_ratio_info,stim_color,negative_vel_threshold,ellipse_ratio_threshold,folder_index, track_index)
%%%% This function plots the velocity and ellipse ratio
    
    light_pulse = area([(stim_while_turning_peaks) (stim_while_turning_peaks+current_rails_dur)], [0.4 0.4],'facecolor',stim_color,'LineStyle','none');
    alpha(0.1)
    hold on;
    yyaxis left
    plot([(stim_while_turning_peaks-500):1:(stim_while_turning_peaks+current_rails_dur)],dummy_array_with_velocity_info','-r','LineWidth',2,'DisplayName','Velocity')
%     yline_vel_threshold=yline(negative_vel_threshold,'--m','LineWidth',2,'DisplayName','Velocity Threshold');
%     ylim([-.2, .2])
    ax=gca;
    ax.YColor='red';
    ylabel('Velocity')
    hold on;
    yyaxis right
    plot([(stim_while_turning_peaks-500):1:(stim_while_turning_peaks+current_rails_dur)],dummy_array_with_ellipse_ratio_info','-b','LineWidth',2,'DisplayName','Ellipse Ratio')
%     yline_er_threshold=yline(ellipse_ratio_threshold,'--c','LineWidth',2,'DisplayName','Ellipse Ratio Threshold');
%     ylim([2, 8])
    ax=gca;
    ax.YColor='blue';    
    ylabel('Ellipse Ratio')
    hold on;
%     set(get(get(light_pulse,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    %legend('show','Location','best','FontSize',12,'NumColumns',2);
    xlim([stim_while_turning_peaks-500 stim_while_turning_peaks+current_rails_dur])
    xlabel('Frames')
    ax = gca;
    ax.FontSize = 13;
    title(strcat('folder id= ',num2str(folder_index), {' '},{' '},{' '},...
        'worm id= ', num2str(track_index),{' '},{' '},{' '},...
    'frame id= ', num2str(stim_while_turning_peaks)),'FontWeight','Normal')
    hold off

end