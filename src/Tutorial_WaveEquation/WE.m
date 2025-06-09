%-------------------------------------------------------------------------
% Math PhD 2025: Physical Modeling
% The Wave Equation with Loss and Forcing
% Dr Michele Ducceschi
% University of Bologna
% 09-06-2025
%------------------------------------------------------------------------

clear all
close all
clc

%---------------------------------------------------
% custom parameters

fs  = 44100 ;        %-- sample rate [Hz]
T   = 5 ;            %-- total sim time [s]
c   = 200 ;           %-- wave speed [m/s]
L   = 0.67 ;         %-- string length [m]

%- forcing parameters
tw   = 3e-3 ;        %-- total contact duration
F0   = 0.5e2 ;        %-- max forcing amp
xp   = 0.8325 ;        %-- input location [frac]
b    = 0.3 ;           %-- noise mod paramter

%-- output parameters
xL  = 0.34 ;         %-- left channel   [frac]
xR  = 0.57 ;         %-- right channel  [frac]

%-- loss parameters
sig0  = 0.3 ;
sig1  = 3e-4 ;

%-- scheme select
ExpFlag = 0 ; %-- 1 = explitic scheme; else implicit

%-- animation parameters
AnimON      = 1 ;
RefreshRate = 4 ;
%---------------------------------------------------

%---------------------------------------------------
% derived parameters

k           = 1 / fs ;
Ts          = floor(T*fs) ;
tv          = (0:Ts-1)*k ;
fv          = (0:Ts-1)*fs/Ts ;

dx           = c*k ;            %-- grid spacing
if ExpFlag == 1
    dx = sqrt(c^2*k^2+4*sig1*k) ;
end

M           = floor(L/dx) ;     %-- grid subintervals
dx           = L/M ;            %-- adjust grid spacing

%-- forcing
tws = floor(tw*fs) ;

% raised cosine
rc  = zeros(1,Ts) ;
rc(1:tws) = 0.5*F0*(1 - cos(2*pi*(0:tws-1)/tws)) ;
% modulated forcing
noise = rand(1,Ts) ;
fin   = rc.*(1+b*noise) ;
Jin   = spreadinterp(xp,M,dx,1,1) ;


%-- readout
Jout      = zeros(2,M-1) ;
Jout(1,:) = spreadinterp(xL,M,dx,1,2) ;
Jout(2,:) = spreadinterp(xR,M,dx,1,2) ;

% matrix init
BCs       = 1 ;   %-- bc flag for Laplacian Build
D2        = laplacian_build(M,L,BCs) ;
Bp        = (1+sig0*k)*speye(M-1) - sig1*k*D2 ;
B0        = 2*speye(M-1) + c^2*k^2*D2 ;
Bm        = (-1+sig0*k)*speye(M-1) - sig1*k*D2 ;

if ExpFlag == 1
    B0        = 2*speye(M-1) + c^2*k^2*D2 + 2*sig1*k*D2 ;
    Bm        = (-1+sig0*k)*speye(M-1) - 2*sig1*k*D2 ;
end


% state init
ym = zeros(M-1,1) ; y0 = zeros(M-1,1) ;

% out init
out = zeros(2,Ts) ;
%---------------------------------------------------

%---------------------------------------------------
% main loop
tic
for n = 1 : Ts

    if ExpFlag == 1
        yp = (B0*y0 + Bm*ym + k^2*Jin*fin(n)) ./ (1+sig0*k) ;
    else
        yp = Bp \ (B0*y0 + Bm*ym + k^2*Jin*fin(n)) ;
    end


    out(:,n) = Jout * y0 ;

    if AnimON == 1
        if mod(n,RefreshRate) == 0
            plot((1:M-1)*dx,y0,'k') ;
            ylim([-3e-3,3e-3]) ; xlim([0,L]) ;
            tt = sprintf('Elapsed time = %0.2g s', (n-1)*k) ;
            title(tt) ;
            xlabel('string length') ; ylabel('y(x)') ;
            drawnow ;
        end
    end

    ym = y0 ; y0 = yp ;

end
elapsed_time = toc

%-- plot results
close all
subplot(2,1,1)
plot(tv,out(1,:)/1e-3,'k') ; hold on ;
plot(tv,out(2,:)/1e-3,'g') ;
legend('left ch', 'right ch') ;
xlabel('$t$ (s)','interpreter','latex') ;
ylabel('$y(t)$ (mm)','interpreter','latex') ;
set(gca,'ticklabelinterpreter','latex','fontsize',16) ;

subplot(2,1,2)
spcL   = abs(fft(out(1,:))) ;
spcR   = abs(fft(out(2,:))) ;
plot(fv,20*log10(spcL),'k') ; hold on ;
plot(fv,20*log10(spcR),'g') ;
legend('left ch', 'right ch') ;
xlabel('$f$ (Hz)','interpreter','latex') ;
ylabel('$\hat y$ (dB)','interpreter','latex') ;
set(gca,'ticklabelinterpreter','latex','fontsize',16) ;
xlim([0,5000]) ;


%-- play sound
vel = diff(out) ;
vel = vel/max(max(abs(vel))) ;
soundsc(vel,fs) ;










