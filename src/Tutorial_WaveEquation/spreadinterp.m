function J = spreadinterp(x_p,M,h,order,type)


%-- x_p  = input point (as a fraction of L), x_p \in [0,1]
%--   M  = number of grid subintervals
%--   h  = grid spacing
%--order = order of lagrange interpolator [0,1,2]
%-- type = 1: spreading 2: interpolator


mp      = floor(M*x_p) ;
alpha   = M*x_p - mp ;
J       = zeros(M-1,1) ;

if order == 0 
    J(mp) = 1 ;
elseif order == 1 
    J(mp-1) = 0.5*(1-alpha) ;
    J(mp+1) = 0.5*(1+alpha) ;
elseif order == 2
     J(mp-1) = 0.5*(alpha-1)*alpha ;
     J(mp)   = (1-alpha)*(1+alpha) ;
     J(mp+1) = 0.5*(alpha+1)*alpha ;
end

if type == 1
    J = J/h ; 
elseif type == 2
    J = J' ;
end
