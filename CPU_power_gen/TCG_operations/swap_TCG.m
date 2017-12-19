function [ Gh, Gv ] = swap_TCG( idx1, idx2, Gh, Gv )

    temp=Gh(idx1,:);
    Gh(idx1,:)=Gh(idx2,:);
    Gh(idx2,:)=temp;
    
    temp=Gh(:,idx1);
    Gh(:,idx1)=Gh(:,idx2);
    Gh(:,idx2)=temp;
    
    temp=Gv(idx1,:);
    Gv(idx1,:)=Gv(idx2,:);
    Gv(idx2,:)=temp;
    
    temp=Gv(:,idx1);
    Gv(:,idx1)=Gv(:,idx2);
    Gv(:,idx2)=temp;
 
end

