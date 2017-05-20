function [ Signal_Rx,Signal_Tx,gain ] = Channel_HW( Signal_Tx,Frec,DAC_port,DAC_FS,Fs,ADC_Device,ADC_port,ADC_FS,T_ADC_iteration_sec )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here: The function operates NI DAQ devices;
%   A/D and D/A in the following order:
% 1) opens sessions with A/D and D/A and defines their parameters (Fs,
% dynamic ranges etc.)
% 2) defines a listener object to the event of the end of the reading of
% the data from the A/D; it triggers a function (Store_data) that takes the
% data read from A/D and saves it to a global variable (BBB) that is
% , of course, available for the entire code 
% and is passed on to the Rx chain as "Signal Rx" (bottom of function)

global BBB


gain=(DAC_FS/max(abs(Signal_Tx)));
Signal_Tx=Signal_Tx*gain; % D/A's output range is [-10,10] volt. we therefore increase output voltage to reach [-9,9] range


%% NI configuration

%%% ADC
aiSession = daq.createSession('ni');
aiSession.Rate = Fs;

if ~isempty(T_ADC_iteration_sec)
    aiSession.DurationInSeconds = T_ADC_iteration_sec;
else
    aiSession.DurationInSeconds = 10;
end

switch ADC_Device
    case '6212'

        Channel_ADC=aiSession.addAnalogInputChannel('Dev1',ADC_port,'Voltage');
    case '9215'
        
        ADC_port=strcat('ai',num2str(ADC_port));
        Channel_ADC=aiSession.addAnalogInputChannel('cDAQ1Mod1',ADC_port,'Voltage');
end

%Channel_ADC=aiSession.addAnalogInputChannel('cDAQ1Mod1','ai0','Voltage');

Channel_ADC.Range=[-ADC_FS,ADC_FS];
lh = addlistener(aiSession,'DataAvailable', @Store_Data);

%%% DAC
aoSession = daq.createSession('ni');
aoSession.Rate = Frec;

Channel_DAC=aoSession.addAnalogOutputChannel('Dev1', DAC_port, 'Voltage');


%% Tx& Rx
queueOutputData(aoSession,Signal_Tx);
%prepare(aiSession)
startBackground(aiSession); %% ADC

%tic
startForeground(aoSession); %% DAC
%toc

%tic
wait(aiSession)
%toc

Signal_Rx=BBB; % global variable modified by the Store_Data function


%% Releases
delete(lh)

end

