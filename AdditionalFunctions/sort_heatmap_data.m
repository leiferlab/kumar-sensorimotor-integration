function [sorted_value,sorted_idx] = sort_heatmap_data(data_to_be_plotted,begin_index,end_index,parameter,sort_order,comparision_method)   
% [sorted_value,sorted_idx] = sort_heatmap_data(data_to_be_plotted,begin_index,end_index,parameter,sort_order,comparision_method)   
% This function will sort the imagesc data. data_to_be_plotted is the raw matrix, 
% begin index, end index, are the two indices between which the function will sort data
% paramter is what you want to sort the data on. For e.g. maybe you want to sort on reversals (i.e. 3)
% sort_order is 'ascend' or 'descend'. Finally comparision method takes 1, 2, or 3 as input.
% 1 is <=, 2 is ==, and 3 is >=
    
    index=[];
    
    for i=1:size(data_to_be_plotted,1)
        if comparision_method==1
            dummy=find(data_to_be_plotted(i,begin_index:end_index)<=parameter);
        elseif comparision_method==2
            dummy=find(data_to_be_plotted(i,begin_index:end_index)==parameter);
        else
            dummy=find(data_to_be_plotted(i,begin_index:end_index)>=parameter);
        end
        if isempty(dummy)
            index(i)=0;
            continue
        end
        index(i)=dummy(1);
    end
    
    [sorted_value,dummy_index]=sort(index,sort_order);
    if isempty(find(sorted_value,1))
        sorted_idx=[1:1:size(data_to_be_plotted,1)];
    else
        sorted_idx=[dummy_index(find(sorted_value,1):end) dummy_index(1:(find(sorted_value,1)-1))];
    end
end