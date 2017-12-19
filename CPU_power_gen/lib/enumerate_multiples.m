function [ multiples ] = enumerate_multiples( x )

    multiples=[];

    for i=1:x
        if mod(x,i)==0
            multiples=[multiples i];
        end
    end

end