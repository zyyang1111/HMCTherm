function [ Gh, Gv, w, h ] = random_swap_TCG( Gh, Gv, w, h )

    n=length(w);
    
    idx=randperm(n,2);
    idx1=idx(1);
    idx2=idx(2);
    
    [ Gh, Gv ] = swap_TCG( idx1, idx2, Gh, Gv );
    
%     display(['swap ' num2str(idx1) ' and ' num2str(idx2)]);

end

