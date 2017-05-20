function [ Signal_Rx,Signal_Tx,gain ] = Channel_HW_w_synch( Signal_Tx,Frec,DAC_port,DAC_FS,Fs,ADC_Device,ADC_port,ADC_FS,F_if,BW,T_ADC_iteration_sec )
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

clear global BBB
clear global CCC
clear global peaks_actual
clear timeout


global BBB
global CCC
global peaks_actual
global timeout

timeout=100; % inside synchro. in case stuck inside

switch BW
    case 5e3
        warning('BW=5kHz does not synchronize well')
        
        if Frec==100e3;
            
            prominence_ref=1.6e4;
            Pref=1.3e9;
            
            %         elseif Frec==250e3
            %             prominence=6.4e3;
            %             Pref=3.64e8;
        end
        
        Width=3;
        BW_chirp=5e3;
        
    case 10e3
        if Frec==100e3;
            
            prominence_ref=1.2e4;
            Pref=7.1e8;
            
        elseif Frec==250e3
            
            prominence_ref=5e3;
            Pref=1.5e8;
        end
        
        Width=3;
        BW_chirp=10e3;
        
    case 20e3 % the chirp signal height decreases on high frequencies somehow, and harms the autocorrelation. Thus, we preferred generating a 10kHz chirp instead of 20kHz
        if Frec==100e3;
            
            prominence_ref=1.2e4;
            Pref=7.1e8;
            
            %             prominence_ref=1.0e4;
            %             Pref=2.7e8;
            %             Width=2.65;
            
        elseif Frec==240e3 % needs to be a multiple of BW
            error('I was unable to find the right parameters')
            
            prominence_ref=2.5e3;
            Pref=1.5e8;
            Width=2.2;
        end
        
        BW_chirp=10e3;
        Width=3;
        
        %         error('prominence was calculated only to 5kHz or to 10kHz BW and to 100kHz Fs and 250kHz or 100kHz Frec')
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
h = fliplr(x); % matched= flipped
x=x*max(abs(Signal_Tx));


if Fs~=Frec
    [p,q]=rat(Fs/Frec);
    h=resample(h,p,q); % chirp Tx was generated with Frc. resampling the filter to be matched to the signal sampled at Fs
    
    % LPF used for fiteribg out the aliased D/A's replica in case Frec=250kHz
    Fpass=(F_if+BW/2)*1.1; % actually, it is F_chip/2.  After upsampling, the right edge of the first replica decreases by N_upsample
    Fstop= Fpass*1.2;
    Hpost_SPEC=fdesign.lowpass('N,Fp,Fst,Ap',100,Fpass,Fstop,0.02,Fs);% in the case of Nguard bad right=10
    Hpost_flt=design(Hpost_SPEC,'equiripple');
else
    Hpost_flt=dfilt.dffir(1);
end
%fvtool(Hpost_flt)


%% Concatenation to the Signal
Signal_Detection_start=[x';zeros(2000,1)];
Signal_Detection_end=[zeros(2000,1);x';zeros(2000,1)];
Signal_Tx=[Signal_Detection_start; Signal_Tx; Signal_Detection_end];

%% Amplification to DAC's Full Scale

gain=(DAC_FS-1)/max(abs(Signal_Tx));
%%% temp
gain=gain/4;
%%%
Signal_Tx=Signal_Tx*gain; % D/A's output range is [-10,10] volt. we therefore increase output voltage to reach [-9,9] range


%% NI configuration

%%% ADC
aiSession = daq.createSession('ni');
aiSession.Rate = Fs;

if ~isempty(T_ADC_iteration_sec)
    aiSession.DurationInSeconds = T_ADC_iteration_sec;
    timeout_background=2*T_ADC_iteration_sec;
    lh = addlistener(aiSession,'DataAvailable', @Store_Data);
else
    aiSession.IsContinuous=true;
    timeout_background=20;
    % lh = addlistener(aiSession,'DataAvailable', @(src,event)Synchro(src,event,prominence,Hpost_flt,h));
    lh = addlistener(aiSession,'DataAvailable', @(src,event)Synchro(src,event,prominence_ref,Hpost_flt,Pref,Width,h));
    
end

switch ADC_Device
    case '6212'
        
        Channel_ADC=aiSession.addAnalogInputChannel('Dev1',ADC_port,'Voltage');
    case '9215'
        
%        ADC_port=strcat('ai',num2str(ADC_port));
         Channel_ADC_0=aiSession.addAnalogInputChannel('cDAQ1Mod1','ai0','Voltage');
         Channel_ADC_1=aiSession.addAnalogInputChannel('cDAQ1Mod1','ai1','Voltage');
         Channel_ADC_2=aiSession.addAnalogInputChannel('cDAQ1Mod1','ai2','Voltage');
         Channel_ADC_3=aiSession.addAnalogInputChannel('cDAQ1Mod1','ai3','Voltage');
        %  aiSession.addAnalogInputChannel('cDAQ1Mod1',ADC_port,'Voltage');
        
end


Channel_ADC_0.Range=[-ADC_FS,ADC_FS];
Channel_ADC_1.Range=[-ADC_FS,ADC_FS];
Channel_ADC_2.Range=[-ADC_FS,ADC_FS];
Channel_ADC_3.Range=[-ADC_FS,ADC_FS];

%%% DAC
aoSession = daq.createSession('ni');
aoSession.Rate = Frec;

Channel_DAC=aoSession.addAnalogOutputChannel('Dev1', DAC_port, 'Voltage');


%% Tx& Rx

queueOutputData(aoSession,Signal_Tx);
%aiSession.NotifyWhenDataAvailableExceeds=8*length(x);
aiSession.NotifyWhenDataAvailableExceeds=6*length(x); % very important; the chirp sequence might be split between 2 consecutive A/D frames (by the size of NotifyWhenDataAvailableExceeds
%)in order to
%reduce the odds of
%such a split we
%set the frame size
%to at least 20
%times the chirp
%size
prepare(aoSession);
prepare(aiSession);

startBackground(aiSession); %% ADC


tic
display('Transmission started')
startForeground(aoSession); %% DAC
wait(aiSession,timeout_background)



Signal_Rx=BBB; % global variable modified by the Store_Data function


%% Releases
aiSession.stop
aiSession.release

aoSession.release

delete(lh)

% plot(CCC)
% [pks,locs,widths,proms] =findpeaks(CCC,'MinPeakProminence',prominence,'MinPeakWidth',3,'Annotate','extents');
%

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
    
    figure(2)
    set(gcf,'windowstyle','docked')
    plot(conv(h,Signal_Tx))
    hold on
    plot(conv(h,Signal_Rx))
    legend(['Tx'],['Rx'])
    
    figure(3)
    set(gcf,'windowstyle','docked')
    findpeaks(conv(h,Signal_Tx),'MinPeakProminence',prominence_ref,'Annotate','extents');
    figure(4)
    set(gcf,'windowstyle','docked')
    findpeaks(conv(h,Signal_Rx),'MinPeakProminence',prominence_ref,'Annotate','extents');
    
    % prominence=prominence*(max(conv(h,Signal_Rx))/3.66e3)
    
    [pks,locs,widths,proms] =findpeaks(conv(h,Signal_Rx),'MinPeakProminence',prominence_ref,'MinPeakWidth',Width,'Annotate','extents')
    
    
    
    %[pks,locs,widths,proms] =findpeaks(conv(h,Signal_Rx),'MinPeakProminence',prominence_ref,'MinPeakHeight',Width,'Annotate','extents')
    
    %     pks
    %     proms
    %     widths
    
    figure(5)
    set(gcf,'windowstyle','docked')
    plot(Signal_Rx)
    title('Acquired Signal')
    
end

%% Chopping the suffix chirp

Signal_Rx_Ref=Signal_Rx(:,1); % Signal_Rx(:,1) is the synchronizing channel
y=conv(h,Signal_Rx_Ref);

%y=conv(h,Signal_Rx);
% locs(1)=max(find(y(1:round(end/2))>peaks_actual(1)*0.95));
% locs(2)=max(find(y(round(end/2)+1:end)>peaks_actual(2)*0.95)+round(length(y)/2));

locs(1)=find(y(1:round(end/2))>peaks_actual(1)*0.95,1); % location of the 1st peak
locs(2)=find(y(round(end/2)+1:end)>peaks_actual(2)*0.95,1)+round(length(y)/2); % location of the 2nd peak

Signal_Rx=Signal_Rx(1:locs(2)-3200,:); % chopping the suffix chirp- located 3200 samples approximately behind the 2nd (last) chirp's autocorrelation peak


end

