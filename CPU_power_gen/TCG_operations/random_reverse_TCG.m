function [ Gh, Gv, w, h ] = random_reverse_TCG( Gh, Gv, w, h )

    n=length(w);
    
    red_idx=[];
    while isempty(red_idx)
        idx=randi(n,1);
        sel_h0_v1=randi(2,1)-1;

        [ Gh_red_idx, Gv_red_idx ] = find_reduction_edges( idx, Gh, Gv );

        if sel_h0_v1==0
            red_idx=Gh_red_idx;
            sel_str='Gh';
        elseif sel_h0_v1==1
            red_idx=Gv_red_idx;
            sel_str='Gv';
        else
            display('ERROR: incorrect value of sel_h0_v1')
        end
    end

    idx1=idx;
    idx2=red_idx(randperm(length(red_idx),1));
    
    [ Gh, Gv ] = reverse_TCG( idx1, idx2, sel_h0_v1, Gh, Gv );
    
%     display(['reverse ' num2str(idx1) ' and ' num2str(idx2) ' in ' sel_str]);
    
end

