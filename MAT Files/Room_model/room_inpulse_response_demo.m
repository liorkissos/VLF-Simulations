clear all; close all; clc;

[sig, Fs]=audioread('2-11 The Voice Of Peace.wma');
speaker=sig(:,1)/max(abs(sig(:,1)));
t=0:1/Fs:(length(speaker)-1)/Fs;
figure
subplot(3,1,1);
plot(t,speaker);
title('Loudspeaker signal');
p = audioplayer(speaker,Fs);
play(p);
fprintf('press enter when playback is done');
%pause
pause(1)

% ---- Load channel model ---------
load('room_impulse_response.mat');
% parameters: 
% Sound velocity 340[m/s]
% Microphone position [1.5 0.5 1] [x y z] [m]
% Louspeaker position [1.5 4.5 1] [x y z] [m]
% Room dimensions [3 5 2.8]       [x y z] [m]
% Reverberation time 0.4 [s]
% Number of taps in room model 8192;
subplot(3,1,2);
plot(room_impulse_response,'g')
title('Room Impulse Response');
% ----------------------------------

% --- Pass the signal through the channel model ------
mic=filter(room_impulse_response, 1, speaker);
% ----------------------------------------------------


subplot(3,1,3);
plot(mic,'r');
title('Microphone signal');

mic=mic/max(abs(mic));
p = audioplayer(mic,Fs);
play(p);