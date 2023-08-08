function test_results=display_significance_results(p_value)
%%%% displays the results of z-test
    
    if isnan(p_value)
        test_results='';
    end
    
    if p_value<0.05 && p_value>=0.01
        test_results='*';
    elseif p_value<0.01 && p_value>=0.001
        test_results='**';
    elseif p_value<0.001 
        test_results='***';
    elseif p_value>0.05 
        test_results='n.s.';
    end

end