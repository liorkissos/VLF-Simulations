function gain=DAC_HW( Signal_Tx,Frec,DAC_port,DAC_FS,Vpeak_out,F_if,BW)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

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


%debug_flag=1;
debug_flag=0;



switch BW
    case 5e3

BW_chirp=5e3;
        
    case {10e3,20e3}

BW_chirp=10e3;
        
    otherwise
        error('prominence was calculated only to 5kHz or to 10kHz BW and to 100kHz Fs and 250kHz or 100kHz Frec')
end
        
Trec=1/Frec;



%% Chirp generation
 
% parameters
t1=2000*Trec;
f0=F_if-BW_chirp/2;
f1=F_if+BW_chirp/2;

% Signal
t=0:Trec:(round(t1/Trec)-1)*Trec;
x=chirp(t,f0,t1,f1,'linear');

% Matched filter
x=x*max(abs(Signal_Tx));


%% Concatenation to the Signal
%Signal_Detection_start=[x';zeros(2000,1)];
Signal_Detection_start=[zeros(2000,1);x';zeros(2000,1)];
Signal_Detection_end=[zeros(2000,1);x';zeros(2000,1)];
Signal_Tx=[Signal_Detection_start; Signal_Tx; Signal_Detection_end];

%% Amplification to DAC's Full Scale

if isempty(Vpeak_out)
    gain=(DAC_FS/max(abs(Signal_Tx)));
else
    gain=(Vpeak_out/max(abs(Signal_Tx)));
end

Signal_Tx=Signal_Tx*gain; % D/A's output range is [-10,10] volt. we therefore increase output voltage to reach [-9,9] range


%% NI configuration

aoSession = daq.createSession('ni');
aoSession.Rate = Frec;

Channel_DAC=aoSession.addAnalogOutputChannel('Dev1', DAC_port, 'Voltage');


%% Tx& Rx

queueOutputData(aoSession,Signal_Tx);

tic
display('Transmission started')
startForeground(aoSession); %% DAC
display('Transmission ended')
toc


%% Releases


aoSession.release

if debug_flag
   
    Signal_Rx=filter(Hpost_flt,Signal_Rx); % filtering the aliased replica
    
    
    %f=linspace(-fs/2,fs/2,length(x));
    f=linspace(-Fs/2,Fs/2,length(x));
    A=fftshift(fft(x));
    figure(1)
    set(gcf,'windowstyle','docked')
    plot(f/1e3,db(abs(A)))
    grid on;grid minor
    xlabel('freq [kHz]')
    title('chirp Tx signal spectrum')
    

end

end

