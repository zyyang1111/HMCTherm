function [ Gh, Gv, w, h ] = random_rotate_TCG( Gh, Gv, w, h )

    n=length(w);
    idx=randi(n,1);
    
    [ w, h ] = rotate_TCG( idx, w, h );
    
%     display(['rotate ' num2str(idx)]);

end

