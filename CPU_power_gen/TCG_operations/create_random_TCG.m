function [ Gh, Gv ] = create_random_TCG( n )

    Gh=sparse(n,n);
    Gv=sparse(n,n);
    connected=eye(n);
    topo=randperm(n);
    edge_ratio_h=rand(1);

    while ~all(reshape(connected,1,[])) %while not all nodes have a connection in some graph
        %%%THIS WHILE LOOP HAS NO GUARENTEE OF CONVERGENCE
        [i j]=find(connected==0); %find all pairs that are unconnected
        k=randi(length(i),1); %choose a random pair
        if find(topo==i(k)) < find(topo==j(k)) %if i(k) is before j(k) in topo
            a=i(k);
            b=j(k);
        else
            a=j(k);
            b=i(k);
        end
        
        %insert edge a->b into some graph
        choose_hv=rand(1);
        if choose_hv<=edge_ratio_h %then the edge goes into Gh
            Gh(a,b)=1;
            Gh_tran=sparse(n,n);
            for i=1:n
                Gh_tran=Gh_tran|Gh^i; %compute transistive closure
            end
            [i j]=find(Gh~=Gh_tran); %find new edges added
            
            %check if any of the transitive edges already exist in other graph
            if ~any(connected(sub2ind([n n],i,j)))
                Gh=Gh_tran; %update the graph to its tran colsure
                accepted=1;
            else
                Gh(a,b)=0; %dont add edge ab
                accepted=0;
            end
            choice_str='h';
        else %then edge goes into Gv
            Gv(a,b)=1;
            Gv_tran=sparse(n,n);
            for i=1:n
                Gv_tran=Gv_tran|Gv^i; %compute transitive closure
            end
            [i j]=find(Gv~=Gv_tran); %find new edges added
            if ~any(connected(sub2ind([n n],i,j)))
                Gv=Gv_tran; %update the graph to its tran colsure
                accepted=1;
            else
                Gv(a,b)=0; %dont add edge ab
                accepted=0;
            end
            choice_str='v';
        end
        if accepted
            connected(a,b)=1;
            connected(b,a)=1;
            connected(sub2ind([n n],i,j))=1;
            connected(sub2ind([n n],j,i))=1;   

%             display(['inserting edge (' num2str(a) ',' num2str(b) ') into graph ' choice_str])
        end
    end

end

