function [ Gh, Gv ] = reverse_TCG( idx1, idx2, sel_h0_v1, Gh, Gv )

    n=size(Gh,1);

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
    G(idx2,idx1)=1;
    
    fout=find(G(idx1,:));
    fin=find(G(:,idx2))';
    
    for k=[fin idx2]
        for l=[fout idx1]
            if G(k,l)==0
                G(k,l)=1;
                other_G(k,l)=0;
                other_G(l,k)=0;
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
    
%     Gh
%     Gv

end

