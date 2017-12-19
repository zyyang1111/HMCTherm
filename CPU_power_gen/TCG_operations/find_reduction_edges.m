function [ Gh_red_idx, Gv_red_idx ] = find_reduction_edges( idx, Gh, Gv )

    Gh_red_idx=[];
    Gv_red_idx=[];

%     size(Gh)
%     idx
    
    Gh_fout=find(Gh(idx,:));
    Gv_fout=find(Gv(idx,:));

    for k=Gv_fout
        k_is_red=1;
        for l=Gv_fout
            Gv_fout_l=find(Gv(l,:));
            if ismember(k,Gv_fout_l)
                k_is_red=0;
            end
        end
        if k_is_red
            Gv_red_idx=[Gv_red_idx k];
        end
    end
    
    for k=Gh_fout
        k_is_red=1;
        for l=Gh_fout
            Gh_fout_l=find(Gh(l,:));
            if ismember(k,Gh_fout_l)
                k_is_red=0;
            end
        end
        if k_is_red
            Gh_red_idx=[Gh_red_idx k];
        end
    end

end