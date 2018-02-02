function [ topologies ] = enumerate_network_topologies( x )

    [ multiples ] = enumerate_multiples( x );
    
    topologies=[];
    
    for i=1:length(multiples)
        for j=1:length(multiples)
            for k=1:length(multiples)
                if multiples(i)*multiples(j)*multiples(k)==x
                    topologies=[topologies; multiples(i) multiples(j) multiples(k)];
                end
            end
        end
    end
    
%     topologies(topologies==1)=NaN;
    
    AR=topologies./repmat(min(topologies,[],2),1,size(topologies,2));
    good_choices=find(max(AR,[],2)<=4);
    
    topologies=topologies(good_choices,:);

end

