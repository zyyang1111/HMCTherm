function [ Gh, Gv ] = create_trivial_TCG( n )

Gh=sparse(n,n);
Gv=sparse(tril(ones(n),-1));

end

