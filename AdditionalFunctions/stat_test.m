function p_value=stat_test(data_1,data_2)
    %%%% data_1 and data_2 should be column vector with one column
    [mean_data_1,~,~]=bootstrap_mean_and_ci(10000,0.05,data_1(:,1));
    [mean_data_2,~,~]=bootstrap_mean_and_ci(10000,0.05,data_2(:,1));

    x1_bar_expt=mean_data_1;
    x2_bar_expt=mean_data_2;
    n1=size(data_1,1);
    n2=size(data_2,1);

    p_hat=(n1*x1_bar_expt+n2*x2_bar_expt)/(n1+n2);
    t_value=(x1_bar_expt -x2_bar_expt)/sqrt((p_hat*(1-p_hat)*(1/n1+1/n2)));
    df_value=n1+n2-2;
    p_value=2*(1-tcdf(abs(t_value),df_value));
end