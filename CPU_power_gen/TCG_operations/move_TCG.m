function [ Gh, Gv ] = move_TCG( idx1, idx2, sel_h0_v1, Gh, Gv )

    %choose selected graph to operate on
    if sel_h0_v1==0
        G=Gh;
        other_G=Gv;
        sel_str='Gh';
    elseif sel_h0_v1==1
        G=Gv;
        other_G=Gh;
        sel_str='Gv';
    end
    
    %verify that idx1,idx2 is a reduction edge in G
    [ G_red_idx, other_G_red_idx ] = find_reduction_edges( idx1, G, other_G );
    if ~ismember(idx2, G_red_idx)
        display(['ERROR: (' num2str(idx1) ',' num2str(idx2) ') is not a reduction edge in graph ' sel_str])
        return
    end
    
    G(idx1,idx2)=0;
    other_G(idx1,idx2)=1;
    
    fout=find(other_G(idx2,:));
    fin=find(other_G(:,idx1))';
    
    for k=[fin idx1]
        for l=[fout idx2]
%             [k l]
            if other_G(k,l)==0
                other_G(k,l)=1;
                G(k,l)=0;
                G(l,k)=0;
            end
        end
    end
    
    %write back to selected graph
    if sel_h0_v1==0
        Gh=G;
        Gv=other_G;
    elseif sel_h0_v1==1
        Gv=G;
        Gh=other_G;
    end

end

