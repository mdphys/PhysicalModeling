%--------------------------------------------------------------------------
%                    Physical Modelling
%               Harmonic Oscillator: Convergence
%                  Dr Michele Ducceschi
%                  University of Bologna
%                        12 May 2025
%--------------------------------------------------------------------------



clear all
close all
clc

%---------------------------------------------------------------
% custom parameters
fsvec   = 0.5*[1,5,10,50,100,500]*1e3 ;
om0     = 1000 ;
T       = 1.0 ; %-- time in seconds 
disp0   = 1 ;
vel0    = 1 ;
AccOrd  = 4 ; %-- accuracy order of ICs
%---------------------------------------------------------------

kvec    = 1 ./ fsvec ;
ErrVec  = zeros(length(kvec),1) ;

for m = 1 : length(kvec)
    
    %-----------------------------
    % derived parameters
    k       = kvec(m) ;
    cs      = cos(om0*k) ;
    Ts      = floor(T / k) ;
    
    %-- initial conditions
    um      = disp0  ;
    if AccOrd == 1
        u0      = um + k*vel0 ;
    elseif AccOrd == 2
        u0 = um + k*vel0 - 0.5*k^2*om0^2*um ;
    elseif AccOrd == 3
        u0 = um + (k*vel0 - 0.5*k^2*om0^2*um)/(1 + k^2 / 6 * om0^2) ;
    elseif AccOrd == 4
        u0 = um + (k*vel0 - 0.5*k^2*om0^2*um - k^4/24*om0^4*um)/(1 + k^2 / 6 * om0^2) ;
    end
    %----------------------------------
    % exact initial conditions (sampled from the analytic solution)
    % um = disp0 ;
    % u0 = disp0 * cos(om0*k) + vel0/om0 * sin(om0*k) ;
    
    %----------------------------------
    % init
    out     = zeros(Ts,1) ; out(1) = um ; out(2) = u0 ;
    tv      = (0 : Ts - 1) * k ;
    exacSol = disp0 * cos(om0*(Ts-1)*k) + vel0/om0 * sin(om0*(Ts-1)*k) ;
    
    %+++++++++++++++++++++++++++++++++++
    
    %+++++++++++++++++++++++++++++++++++
    % main loop
    
    for n = 3 : Ts
        
        up      = 2*cs*u0- um     ;
        out(n)  = up ;
        um      = u0 ;
        u0      = up ;
        
    end
    
    ErrVec(m) = exacSol -  up   ;
    
end
%------------------------------------
% plots

loglog(fsvec,abs(ErrVec),'linewidth',3) ; grid on ; grid minor ;
set(gca, 'fontsize',20)
xlabel('fs (Hz)'); ylabel('|Err|')
