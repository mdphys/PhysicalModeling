%---- An additive synthesiser based on the SHO
%     Michele Ducceschi
%     12 May 2025
%-----------------------------------------------------

function [yMono,y] = additivesynt(Om0,T60s,b,fin,Tadd,fs)


N    = length(Om0) ;
Ts   = length(fin) + floor(Tadd*fs) ;
sigs = 3*log(10)./ T60s ;
k    = 1/fs ;
bT   = b' ;

yMono = zeros(Ts,1) ;
y     = zeros(N,Ts) ;

ym    = zeros(N,1) ;
y0    = zeros(N,1) ;

exps  = exp(-sigs*k) ;
P0    = 2*exps.*cos(Om0*k) ;
Pm    = exps.^2 ;
Pf    = k^2*exps ;



for n = 1 : Ts 

    fext = 0 ;
    if n < length(fin)
        fext = fin(n) ;
    end
     
    % EXACT OSCILLATOR 
    yp = P0.*y0 - Pm.*ym + Pf*fext ;
    
    y(:,n) = yp' ;
    yMono(n) = bT * yp ;

    ym = y0 ; y0 = yp ;

end


yMono = yMono / max(abs(yMono)) ;
