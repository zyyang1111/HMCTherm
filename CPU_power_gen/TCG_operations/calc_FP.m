function [ x, y ] = calc_FP( Gh, Gv, w, h )

    n=size(Gh,1);

    %calc indeg and outdeg of each node in each graph
    Gv_outdeg=sum(Gv,2);
    Gv_indeg=sum(Gv,1);

    Gh_outdeg=sum(Gh,2);
    Gh_indeg=sum(Gh,1);

    %scale nodes by -w/-h and find shortest path (i.e. longest path in w/h)
%     Gh
%     Gv
    Gh_dist=all_shortest_paths(Gh.*repmat(-w,1,n));
    Gv_dist=all_shortest_paths(Gv.*repmat(-h,1,n));

    %find the distance from each sorce (indeg=0) node to each other node
    Gh_source_dist=Gh_dist(find(Gh_indeg==0),:);
    Gv_source_dist=Gv_dist(find(Gv_indeg==0),:);

    %x,y are the maximum path length accross all sorce nodes
    %multiply path lengths by -1 to get positive lengths
    x=max(-Gh_source_dist,[],1)';
    y=max(-Gv_source_dist,[],1)';

end