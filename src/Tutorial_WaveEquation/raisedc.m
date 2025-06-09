clear all
close all

%-------------------------------------------
%-- custom parameters
fs  = 44100 ;      %-- sample rate [Hz]
T   = 0.2 ;        %-- duration [s]

% forcing params
F0  = 1 ;          %-- amplitude 
tw  = 1e-3 ;       %-- contact duration [s]
b   = 0.1 ;        %-- noise mod parameter
%-------------------------------------------

%-------------------------------------------
%-- post proc
Ts  = floor(T*fs) ;
tws = floor(tw*fs) ;
tv  = (0:Ts-1)/fs ;

% raised cosine
rc  = zeros(1,Ts) ;
rc(1:tws) = 0.5*F0*(1 - cos(2*pi*(0:tws-1)/tws)) ;

% modulated forcing 
noise = rand(1,Ts) ;
f     = rc.*(1+b*noise) ;

plot(tv,rc,tv,f) ;  xlim([0,1.5e-3]) ;
set(gca,'fontsize',20)
xlabel('t (s)'); ylabel('f (t)'); 
legend('raised cos', 'modulated raised cos')



