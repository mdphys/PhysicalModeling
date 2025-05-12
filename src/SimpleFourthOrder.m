%--------------------------------------------------------------------------
%                    Physical Modelling
%               Five-point time stencil for a particle in a free field
%                        (UNSTABLE!!)
%                  Dr Michele Ducceschi
%                  University of Bologna
%                        12 May 2025
%--------------------------------------------------------------------------

clear all
close all
clc

%--- custom parameters

fs    = 1000 ;    %- sample rate [Hz]
T     = 0.1 ;     %- total sim time [s]
disp0 = 1 ;       %- initial displacement 
vel0  = 1 ;       %- initial velocity
%--------------------------------------------------------------------------




%--- derived parameters
 
k = 1/fs ; %- sample rate

%- ICs 
u0 = disp0 ;
um = u0 - k*vel0 ;
umm = um - k*vel0 ;
up = u0 + k*vel0 ;

Ts = floor(T*fs) ;
Ts = 20 ;

out = zeros(Ts,1) ;  %- init out 
%--------------------------------------------------------------------------


%--- main loop

for n = 1 : Ts 

    %- implement fourth-order time difference for the second derivative

    upp = 16*up - 30*u0 + 16*um - umm ;

    out(n) = u0 ;

    umm = um ; um = u0 ; u0 = up ; up = upp ;

end


plot(out)
  