function [ Gh, Gv ] = create_TCG_from_FP( x, y, w, h )

    n=length(x);
    Gh=sparse(n,n);
    Gv=sparse(n,n);
    connected=eye(n);

    for i=1:n
        for j=1:n
            if x(i)>x(j) %if j is to the left of i
                if y(i)<y(j)+h(j) && y(j)<y(i)+h(i) %i and j overlap on y axis
                    Gh(j,i)=1;
                    connected(j,i)=1;
                    connected(i,j)=1;
                end 
            end

            if y(i)>y(j) %if j is below i
                if x(i)<x(j)+w(j) && x(j)<x(i)+w(i) %i and j overlap on x axis
                    Gv(j,i)=1;
                    connected(j,i)=1;
                    connected(i,j)=1;
                end
            end
        end
    end
    
%     connected

    %%%transitively close graphs
    Gh_tran=sparse(n,n);
    Gv_tran=sparse(n,n);

    for i=1:n
        Gh_tran=Gh_tran|Gh^i;
        Gv_tran=Gv_tran|Gv^i;
    end
    
    [i j]=find(Gh~=Gh_tran); %find new edges added
    connected(sub2ind([n n],i,j))=1;
    connected(sub2ind([n n],j,i))=1;

    [i j]=find(Gv~=Gv_tran); %find new edges added
    connected(sub2ind([n n],i,j))=1;
    connected(sub2ind([n n],j,i))=1;
    
    Gh=Gh_tran;
    Gv=Gv_tran;

    %any unconnected nodes should be connected in Gh from left to right
    
    [i j]=find(~connected);
    for k=1:length(i)
        if x(i(k))>x(j(k)) %if j is to the left of i
            Gh(j(k),i(k))=1;
            connected(i(k),j(k))=1;
            connected(j(k),i(k))=1;
        end
    end
    
    %transitively close Gh
    %this may or may not be nessessary
    for i=1:n
        Gh=Gh|Gh^i;
    end
end