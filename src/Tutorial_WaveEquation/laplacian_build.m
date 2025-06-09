function D2 = laplacian_build(M,L,bc)

%---- this function creates the 1D Laplace operator over a domain of length
%---- L, using M grid intervals, and using 3 different kinds of boundary
%---- conditions: 1: Dir Dir, 2: Neu Neu, 3: Dir Neu


h = L/M ;
v = ones(M+1,1) ;


if bc == 1
    D2 = spdiags([v,-2*v,v],-1:1,M-1,M-1) ;
elseif bc == 2
    D2 = spdiags([v,-2*v,v],-1:1,M+1,M+1) ;
    D2(1) = -1 ;
    D2(end) = -1 ;
elseif bc == 3
    D2 = spdiags([v,-2*v,v],-1:1,M,M) ;
    D2(end) = -1 ;
end


D2 = D2 / h^2 ;

end


