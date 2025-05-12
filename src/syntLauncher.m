clear all
close all
clc

N           = 10 ; % number of oscillators 
Om0         = rand(N,1)*10000*2*pi ; % frequencies
T60s        = 0.2 + 10.0*rand(N,1) ; % decay times
b           = -0.5 + rand(N,1) ;
fs          = 48000 ;
Tadd        = 10 ; %- added time in second after forcing has stopped
Tf          = 0.05 ;
fin         = sawtooth((0:floor(Tf*fs))*2*pi*200/fs,0.5) ; % input forcing
fin         = fin'.*hanning(length(fin)) ;
fin = zeros(floor(Tf*fs),1) ; fin(10) = 1 ;
[yMono,y]   = additivesynt(Om0,T60s,b,fin,Tadd,fs) ;

%soundsc(yMono,fs) ;

