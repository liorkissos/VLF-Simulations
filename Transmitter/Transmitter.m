function [Signal_Tx,Fs,OFDM_config,Testing_data,IF_chain_config ] = Transmitter(Symbol_stream_Tx,Coding_config,OFDM_config,IF_chain_config,Simulation_config )
%function [Signal_Tx,Fs,OFDM_config,Testing_data,N_frames_Tx,R1,R2,F_chip,N_upsample ] = Transmitter(Symbol_stream_Tx, OFDM_config,IF_chain_config,N_symbols,P_total,Simulation_config )



%TRANSMITTER Summary of this function goes here
%   Detailed explanation goes here

%%% Dummy # 3: TEST the BANCH BROWSER

%% Config
%test_synch_channel_flag=0; % test only synchronization, impulse response and channel estimation
%test_signal_processing_flag=1;
test_signal_processing_flag=0;

Group_delay_total=0;

%% Input Config

Coding_flag=Coding_config.Coding_flag;
N_coding=Coding_config.N_coding;
K_coding=Coding_config.K_coding;

F_chip =OFDM_config.F_chip;
N_FFT = OFDM_config.N_FFT;
N_data = OFDM_config.N_data;
Npilots   =OFDM_config.Npilots;
Amp_pilots_dB  = OFDM_config.Amp_pilots_dB;
Nguard_band_left  = OFDM_config.Nguard_band_left;
Nguard_band_right  = OFDM_config.Nguard_band_right;
N_CP = OFDM_config.N_CP;
P_total=OFDM_config.P_total;
T=OFDM_config.T;
B_signal_net_Hz=OFDM_config.B_signal_net_Hz;
N_preamble_CE=OFDM_config.N_preamble_CE;
N_preamble_synch=OFDM_config.N_preamble_synch;
Enhancement_prmbl_CE=OFDM_config.Enhancement_prmbl_CE;
Enhancement_prmbl_synch=OFDM_config.Enhancement_prmbl_synch;



F_if=IF_chain_config.F_if;
Frec_req=IF_chain_config.Frec_req;
Pre_emphasis_flag=IF_chain_config.Pre_emphasis;

Input_Type=Simulation_config.Input_Type;
Time_sync_flag=Simulation_config.Time_sync_flag;
BB_IF_flag=Simulation_config.BB_IF_flag;
Configuration=Simulation_config.Configuration;
N_symbols=Simulation_config.N_symbols;
Ref_Cfg_flag=Simulation_config.Ref_Cfg_flag;

%% OFDM params


if(mod(N_data+1,Npilots-1)~=0)
    error('N_data_sub_carriers should be an Integer division of Npilots-1') % to enable uniform scattering of pilot subarriers
end

Amp_pilots=db2pow(Amp_pilots_dB);

if(N_CP>N_FFT)
    error('N_CP should be smaller than N_FFT')
end

if mod(log2(N_FFT),1)~=0 %
    error('N_FFT should be a power of 2')
end

%% Tx Signal Generation %%%%%%%%%%%%%%%%%%

%% Frequency Domain %%

%% 1) Power Calculation


P_data=P_total/(Amp_pilots*Npilots+N_data); % average power of data subcarrier

%% 2) Frame construction and Pilots& DC insertion insertion
%%% Burst initialization: each row represents an OFDM symbol ("super symbol")


N_frames_Tx=length(Symbol_stream_Tx)/N_data;

data_subcarrier_matrix_Tx=reshape(Symbol_stream_Tx,[N_data, N_frames_Tx]);% OFDM frames sprawl along the column dimension. e.g: column # 1 of the matrix is the 1st OFDM frame

data_subcarrier_matrix_Tx_net=data_subcarrier_matrix_Tx; % only pure data symbols/samples/subcarriers. needed for later calculation related noise addition

place_low=(N_FFT/2+1)-(Nguard_band_left+Npilots/2)-1; % the place of the DC carrier
data_subcarrier_matrix_Tx=[data_subcarrier_matrix_Tx(1:place_low,:);zeros(1,N_frames_Tx);data_subcarrier_matrix_Tx(place_low+1:end,:)];

Pilots_matrix_Tx(1,:)=sqrt(Amp_pilots*P_data)*ones(1,N_frames_Tx); % initialization of pilots matrix; needed in receiver

if Npilots
    OFDM_matrix_Tx_f=[];
    chunk_size=(N_data+1)/(Npilots-1);
    nn=1;
    for kk=1:Npilots-1 % insertion of pilot (a row at the amplitude of sqrt(Amp_pilots*P_data) at every "chunk" rows
        
        chunk=data_subcarrier_matrix_Tx(nn:nn+chunk_size-1,:);
        OFDM_matrix_Tx_f=cat(1,OFDM_matrix_Tx_f,cat(1,(sqrt(Amp_pilots*P_data))*ones(1,N_frames_Tx),chunk));  % pilot subacarrier power=Amp_pilots*average(data subcarrier power). need to add power factorization of pilots versus data bins
        
        Pilots_matrix_Tx(kk+1,:)=sqrt(Amp_pilots*P_data)*ones(1,N_frames_Tx);
        
        nn=nn+chunk_size;
    end
    
    OFDM_matrix_Tx_f((N_data+1)+Npilots,:)=(sqrt(Amp_pilots*P_data))*ones(1,N_frames_Tx); % insertion of pilots at the righthand side of the OFDM symbol (frequency domain)
    %OFDM_matrix_Tx_f((N_data-1)+Npilots,:)=(sqrt(Amp_pilots*P_data))*ones(1,N_frames_Tx); % insertion of pilots at the righthand side of the OFDM symbol (frequency domain)
    
    
else % without pilots
    OFDM_matrix_Tx_f=data_subcarrier_matrix_Tx;
    
end

E1=sum(sum((abs(data_subcarrier_matrix_Tx_net)).^2)); % Only payload subcarriers Energy, before concatenatin of DC zeros
E2=sum(sum((abs(OFDM_matrix_Tx_f)).^2)); % Total OFDM block (with pilots and guard bands, but w/o CP) Energy, without the preamble_CE
R1=E1/E2; % frequency domain Enery ratio (<1)

%% 3) Preamble CE (channel estimation) insertion


%%% Sequence taken from: IEEE C802.16a-02/17 document, which contains
%%% particularly low PAPR sequences

switch N_FFT
    case 512
        load('PN_seq_512.mat')
    case 256
        load('PN_seq_256.mat')
    case 64
        % load('PN_seq_64.mat')
        load('PN_seq_64_1.mat')
end

Preamble_CE_stream_Tx=AAA.';

P_preamble_CE=P_total/((N_data)+Npilots); % the subcarriers in this symbol all hav ethe same power
Preamble_CE_stream_Tx=Preamble_CE_stream_Tx*sqrt(P_preamble_CE/1); % scaling the preamble_CE subcarriers to P_preamble_CE power; the average power per premable_CE subcarrier. before scaling, since these are BPSK symbols, their power is (+/-1)^2=1
%%% concatenating the CE ref

%N_preamble_CE=10;
Trainining_Sequence_Tx=repmat(Preamble_CE_stream_Tx,[1 N_preamble_CE ]);
OFDM_matrix_Tx_f=cat(2,Trainining_Sequence_Tx,OFDM_matrix_Tx_f); % the 1st preamble_ce.
% OFDM_matrix_Tx_f=cat(2,Preamble_CE_stream_Tx,OFDM_matrix_Tx_f); % the 2nd preamble_ce . % BFB, p.721, equation (20.28): need to be identical to the previous preamble


N_frames_Tx=N_frames_Tx+N_preamble_CE; % reflects the addition of of the preamble_CE frame

OFDM_Ref_symbol_Tx_f=Preamble_CE_stream_Tx; % not needed in priciple

E3=sum(sum((abs(OFDM_matrix_Tx_f)).^2)); % Total OFDM block (with pilots and guard bands, but w/o CP) Energy, without the preamble_CE
R2=E2/E3; % frequency domain Enery ratio (<1)

%% 4) Guard bands insertion

OFDM_matrix_Tx_f_wo_GB=OFDM_matrix_Tx_f;

OFDM_matrix_Tx_f=[zeros(Nguard_band_left,N_frames_Tx);OFDM_matrix_Tx_f;zeros(Nguard_band_right,N_frames_Tx)];


if test_signal_processing_flag
    if abs(OFDM_matrix_Tx_f(:,3)'*OFDM_matrix_Tx_f(:,3)-P_total)>0.2*P_total % every column of the matrix needs to have a total power of P_total
        error('Power constraint not applied')
    end
    
end


%E4=sum(sum((abs(OFDM_matrix_Tx_f_wo_GB)).^2)); % Only payload subcarriers Energy, before concatenatin of DC zeros
E4=sum(sum((abs(OFDM_matrix_Tx_f)).^2)); % Total OFDM block (with pilots and guard bands, but w/o CP) Energy, without the preamble_CE
R3=E3/E4; % should result in 1, aince guard bands do not add any power

%% Time Domain %%

%% 5) IDFT:
%%% the IDFT mutiplies the 1st carrier with 0*n [rad] (n=time index)
%%% and the last carrier with 2*pi*n [rad]. Since we want the 1st to be
%%% multiplied with -pi*n and the last with +pi*n, we need to apply
%%% ifftshift to the frequency domain OFDM symbol priot to the application
%%% of ifft

OFDM_matrix_Tx_f=ifftshift(OFDM_matrix_Tx_f,1); %  fftshift ove the columns. see comment above
OFDM_matrix_Tx_t=ifft(OFDM_matrix_Tx_f); % the famous OFD IDFT operation. The Matlab function includes the 1/N_FFT factor required by OFDM


if test_signal_processing_flag
    figure(1)
    set(gcf,'windowstyle','docked');
    stem(real(OFDM_matrix_Tx_f(:,5)))
    hold on
    stem(imag(OFDM_matrix_Tx_f(:,5)))
    hold off
    grid minor
    title('OFDM subcarriers (after fftshift)')
    
    figure(2)
    set(gcf,'windowstyle','docked');
    B=fftshift(fft(OFDM_matrix_Tx_t(:,5)));
    plot(db(abs(B)))
    grid on
    
end

E5=size(OFDM_matrix_Tx_t,1)*sum(sum((abs(OFDM_matrix_Tx_t)).^2)); % scaling by the column size ( the ifft operation is applied to the columns of the matrix) to match Parseval theorem for DFT
R4=E4/E5; % sanity check; needs to be equal to 1: Parseval. will not be used in the noise calculation
E5=E5/size(OFDM_matrix_Tx_t,1);% scaling back to facilitate the calculation later. the R4 factor has no effect anyway (since equal to 1)

%% 6) CP insertion
Signal_matrix_Tx=cat(1,OFDM_matrix_Tx_t(end-N_CP+1:end,:),OFDM_matrix_Tx_t);


%%% Power ratios
E6=sum(sum((abs(Signal_matrix_Tx)).^2)); % Total OFDM block (with pilots and guard bands, but w/o CP) Energy, without the preamble_CE
R5=E5/E6; %

%% 7) Preamble CE power enhancment

%%% Only premable_CE is amplified here, as its PAPR preserves its
%%% properties along the Tx chain versus data's PAPR. Whereas,
%%% Preamble_synch a the interpolation stage varies its PAPR by much, thus
%%% we enhance it later.

%if strcmp(Configuration,'Operational')
%%% PAPR testing

P_max_preamble_CE=max(conj(Signal_matrix_Tx(:,1)).*Signal_matrix_Tx(:,1));
P_avg_preamble_CE=mean(conj(Signal_matrix_Tx(:,1)).*Signal_matrix_Tx(:,1));
PAPR_preamble_CE=db(P_max_preamble_CE/P_avg_preamble_CE,'power');

P_max_data=max(max(conj(Signal_matrix_Tx(:,1+N_preamble_CE:end)).*Signal_matrix_Tx(:,1+N_preamble_CE:end)));
P_avg_data=mean(mean(conj(Signal_matrix_Tx(:,1+N_preamble_CE:end)).*Signal_matrix_Tx(:,1+N_preamble_CE:end)));
PAPR_data=db(P_max_data/P_avg_data,'power');
%%% power enhancement

if db(P_max_data,'power')-db(P_max_preamble_CE,'power')<Enhancement_prmbl_CE
    error('impossible to enhance preamble CE')
end



Signal_matrix_Tx(:,1:N_preamble_CE)=Signal_matrix_Tx(:,1:N_preamble_CE)*db2mag(Enhancement_prmbl_CE); %enhancement of the preamble_CE by 4 dB. checked and saw that the proposed 256 long sequence has a PAPR lower than the data sequence's always by higher than 4dB
%Signal_matrix_Tx(:,1:N_preamble_CE)=Signal_matrix_Tx(:,1:N_preamble_CE)*db2mag(3.5); %enhancement of the preamble_CE by 4 dB. checked and saw that the proposed 256 long sequence has a PAPR lower than the data sequence's always by higher than 4dB


if test_signal_processing_flag
    %P_preamble_sync=conj(OFDM_preamble_stream_Tx_t).*OFDM_preamble_stream_Tx_t;
    P_preamble_CE=conj(Signal_matrix_Tx(:,2)).*Signal_matrix_Tx(:,2);
    P_data=conj(Signal_matrix_Tx(:,3+N_preamble_CE)).*Signal_matrix_Tx(:,3+N_preamble_CE);
    
    %plot(db(P_preamble_sync,'power')); hold on;
    plot(db(P_preamble_CE,'power'));hold on;plot(db(P_data,'power'),'g')
    title('Instatneous power of various parts of the packet AFTER(!) amplification')
    legend('Preamble synch','Preamble CE','Data (varies at every run)')
    ylabel('dB')
    xlabel('time [Tchips]')
    grid on;grid minor
    set(gcf,'windowstyle','docked')
    
end
% %%% power ratio calculations
E7=sum(sum((abs(Signal_matrix_Tx)).^2)); % Total OFDM block (with pilots and guard bands, but w/o CP) Energy, without the preamble_CE
R6=E6/E7; %

%% 7*) Impulse insertion
%if test_signal_processing_flag

if strcmp(Input_Type,'Impulse')
    
    %%% Creation of Impulses at each OFDM symbol, apart from 2 preamble
    %         %%% symbols
    %         OFDM_matrix_Tx_t=zeros(size(OFDM_matrix_Tx_t));
    %         OFDM_matrix_Tx_t(1,:)=ones(1,size(OFDM_matrix_Tx_t,2)); % The impulse at the head of each OFDM symbol
    %         OFDM_matrix_Tx_t(:,1)=OFDM_preamble_Tx_t;
    %         OFDM_matrix_Tx_t(:,2)=OFDM_preamble_Tx_t; % 2 premble symbols are required for a correct Minn& Zeng sychronization
    %         OFDM_matrix_Tx_t(:,3)=zeros(N_FFT,1); % 2 premble symbols are required for a correct Minn& Zeng sychronization
    
    %%%% 2 preamble symbols, 1 void symbol, 2 impulse symbols then void
    OFDM_matrix_Tx_t=zeros(size(OFDM_matrix_Tx_t));
    OFDM_matrix_Tx_t(:,1)=OFDM_preamble_Tx_t;
    OFDM_matrix_Tx_t(:,2)=zeros(N_FFT,1); % 2 premble symbols are required for a correct Minn& Zeng sychronization
    OFDM_matrix_Tx_t(1,3:4)=P_total*ones(1,1); % The impulse at the head of each OFDM symbol
    
end

%end

%% Tx DSP FE %%%

%% 8) Parallel to Serial

Signal_Tx=reshape(Signal_matrix_Tx,[(N_FFT+N_CP)*(N_frames_Tx),1]); % N_frames_Tx+1: N_frames_Tx data frames and 1 preamble frame

%% 9) Preamble Synch insertion
%%% time domain directly

% N=(2^m)-1 length PN series, constructed by a prime
% polynomial (CCF architecture [canonic representation]), and thereofre is
% a m-length sequence with exceptional Autotcorrelation properties
% see: www.gaussianwaves.com/2010/09/maximum-length-sequences-m-sequences-2
% good correlation properties are needed: 1) speed of convergence of
% algorithm (BFB p. 147, p.155. BLM p. 437) 2) steady state (execess error) value (BFB., p.153)
% it conserves the properties of orthogonality under BPSK modulation


%%% generating the basic subframe

%m=log2(N_FFT/4); % Minn & Zeng. L=N_FFT/4


%N_preamble_synch=8;
m=log2(N_preamble_synch*N_FFT/4); % Minn & Zeng. L=N_FFT/4

switch m % taken from; http://www.mathworks.com/help/comm/ref/pnsequencegenerator.html
    case 4
        initial_conditions = randi([0 1],m,1);
        Polynomial=[4 3 0];
    case 5
        Polynomial=[5 3 0];
    case 6
        Polynomial=[6 5 0];
        load('PN_256_init'); % loading initial_conditions chosen best for time synch
    case 7
        initial_conditions = randi([0 1],m,1);
        Polynomial=[7 6 0];
    case 8
        initial_conditions = randi([0 1],m,1);
        Polynomial=[8 6 5 4 0];
    case 9
        initial_conditions = randi([0 1],m,1);
        Polynomial=[9 5 0];
    case 10
        initial_conditions = randi([0 1],m,1);
        Polynomial=[10 7 0];
        
    otherwise
        error('no polynomial matches required premabole length')
end

N_preamble_symbols=(2^m-1)/log2(2); % preamble is always BPSK, thus log2(2)

hpn = comm.PNSequence('Polynomial',Polynomial,'SamplesPerFrame', N_preamble_symbols,'InitialConditionsSource','Input port');

preamble_Tx = step(hpn,initial_conditions);

hMod_preamble=comm.BPSKModulator;
Preamble_stream_Tx=step(hMod_preamble,preamble_Tx);

if test_signal_processing_flag % test the orthogonality ( based on Gaussinawave book. simulation is in Experiments folder)
    %a=Preamble_stream_Tx(2:end);
    a=Preamble_stream_Tx(1:end);
    Np=2;
    a_present = repmat(a,[Np,1]); %copy the sequence Np times
    L = length(a_present);
    %Compute the autocorrelation
    a_past = a_present'; %delay=0
    for k = 1:L
        C(k) = (a_past*a_present)/L; % Correlation. C(1) should produce maximal autocorrelation
        % cyclic extension of a_past
        a_past_out = a_past(end);
        a_past(2:end) = a_past(1:end-1);
        a_past(1) = a_past_out;
    end
    
    stem(real(C));
    set(gcf,'windowstyle','docked');
    
end

%%% preparing the frame (see Minn& Zheng eq.9)

OFDM_preamble_stream_Tx_t=[Preamble_stream_Tx;Preamble_stream_Tx;-Preamble_stream_Tx;-Preamble_stream_Tx]; % see Minn& Zheng eq.9
OFDM_preamble_stream_Tx_t=[OFDM_preamble_stream_Tx_t;zeros(4,1)]; % since Minn&Zeng preamble always comped of 4 sub sequences, we will always need to add 4 zeros

%%% power constraint application
OFDM_preamble_stream_Tx_t=OFDM_preamble_stream_Tx_t*sqrt(P_total/N_FFT^2)/sqrt(mean(abs(OFDM_preamble_stream_Tx_t).^2));% time domain OFDM symbol's energy is lower by N_FFT (Parseval), and P_total is frequency domain energy, thus P_total/N_FFT is time domain energy.P_total/N_FFT^2 is the energy per subcarrier

%%% concatenatio of a single(!) CP
OFDM_preamble_stream_Tx_t=[OFDM_preamble_stream_Tx_t(end-N_CP+1:end,:);OFDM_preamble_stream_Tx_t];


%%% Concatenation: Preamble_synch to Data

Signal_Tx=[OFDM_preamble_stream_Tx_t;Signal_Tx];

%%% Power ratios
%%% Power ratios
E8=sum((abs(Signal_Tx).^2)); % Total OFDM block (with pilots and guard bands, but w/o CP) Energy, without the preamble_CE
R7=E7/E8; %

Signal_Tx_OFDM=Signal_Tx;

if test_signal_processing_flag
    close all
    figure(1)
    plot(conj(Signal_Tx).*Signal_Tx)
    set(gcf,'windowstyle','docked')
    
    
    figure(2)
    Fs=F_chip;
    set(gcf,'windowstyle','docked');
    f=linspace(-Fs/2,Fs/2,length(Signal_Tx));
    %A=fftshift(fft(Signal_Tx_upsampled,length(f)));
    C=fftshift(fft(Signal_Tx,length(f)));
    % plot(f/1e3,db(abs(A)))
    plot(f/1e3,db(abs(C)))
    %title(['Signal Tx after upsampling by ',num2str(N_upsample),''])
    xlabel('freq [kHz]')
    grid on
    grid minor
end

%% 10) length ratio Calculation: needed for noise addition

if Coding_flag
    R_lengths=length(Symbol_stream_Tx)*(K_coding/N_coding)/(length(Signal_Tx));
else
    R_lengths=length(Symbol_stream_Tx)/(length(Signal_Tx)); % another scaling factor needed for AWGN addition caculation (see Channel_SW function as well as a summary on the subject)
end

%% 11) Interpolation

%%% before Interpolation

% the nyquist zone boundary at the output of the OFDM modulator.
%it operates there at a rate of F_chip. also equals (N_FFT/2)*1/T; %every bin's BW is 1/T. on the positive frequency side there is a total number of N_FFT bins
Fs=F_chip; % current Fs

%%% Interpolation Filter

if strcmp(Configuration,'Operational')  % in any case other than Calibration, the interpolation factor and thus filter need to be flexible in order to adapt to the variable F_chip and the maximum Fsample cnstrained by the NI 6212
    
    N_upsample=floor(Frec_req/F_chip);
    Fs=Fs*N_upsample; % Fnyquist after upsampling. needed since the filter follows the upsampling. limited to 250ksps due to D/A (NI 6212) constraints
    
    
    if  Ref_Cfg_flag% operationl and reference configuration
        load('Filters_Ref_Cfg.mat')
        
    else
        
        Fpass=B_signal_net_Hz; % actually, it is F_chip/2.  After upsampling, the right edge of the first replica decreases by N_upsample
        Fstop=F_chip/2+min(Nguard_band_left,Nguard_band_right)/T;%since it is a BB signal, the 2nd replica lies at an interval of Nguard_band_left/T from Fnyquist (the current Fs/2).
        
        Fpass=Fpass*1.035; % do not erase! original: Integer group delay
        
        Hinterp_SPEC=fdesign.lowpass('Fp,Fst,Ap,Ast',Fpass,Fstop,0.002,65,Fs);% in the case of nguard band =6
        Hinterp_flt=design(Hinterp_SPEC,'kaiserwin');
        
    end
    
    
elseif strcmp(Configuration,'Calibration')  % only calibration case: use specially designed interpolation filter, adapted to N_upsample=14
    N_upsample=14;
    Fs=min(Frec_req,Fs*N_upsample); % limiting of the Fsample to the maximum possible by the NI 6212, which is 250kHz
    F_chip=F_chip*(F_chip*N_upsample/Fs); % effective F_chip may be different from the requested one due to the limit of Fsample to 250khz
    
    load('Filters.mat') % load all the filters at once
    
    
    if mod(Hinterp_flt.order,2)~=0 && Time_sync_flag==0 % integer delay needed only when no synchronization nor equalization is made
        error('Interpolating filter group delay is not an integer number')
    end
    
end

Signal_Tx=upsample(Signal_Tx,N_upsample);

Signal_Tx_upsampled=Signal_Tx; % for debugging purposes

Signal_Tx=N_upsample*conv(Hinterp_flt.numerator,Signal_Tx,'full');  % conv instead of filter so as not to loose the suffix samples of the signal, that might cause loosing an actual payload frame
Group_delay_interp=Hinterp_flt.order/2;

Group_delay_total=Group_delay_total+Group_delay_interp;

if test_signal_processing_flag
    
    close all
    plot(conj(Signal_Tx).*Signal_Tx)
    set(gcf,'windowstyle','docked')
    
    
    Signal_Tx_interp=Signal_Tx;
    %Signal_Tx_interp1=Signal_Tx1;
    
    Group_delay_post_interpolation=Group_delay_total;
    % Pilots might spread along several bins becuase the bins of FFT applied over the whole
    % signal are narrower than those of the OFDM FFT, and therefor hte energy is spread along several bins
    figure(3)
    set(gcf,'windowstyle','docked');
    f=linspace(-Fs/2,Fs/2,length(Signal_Tx));
    A=fftshift(fft(Signal_Tx_upsampled,length(f)));
    C=fftshift(fft(Signal_Tx,length(f)));
   % plot(f/1e3,db(abs(C)))
    plot(f/1e3,db(abs(A)),f/1e3,db(abs(C)))
    title(['Signal Tx after upsampling by ',num2str(N_upsample),''])
    xlabel('freq [kHz]')
    grid on
    grid minor
    
    A=Signal_Tx_OFDM;
    B=Signal_Tx_interp(Group_delay_post_interpolation+1:N_upsample:end);
    
    %D=Signal_Tx_interp1(1:N_upsample:end);
    %EVM_Lior_Interpolation=db(mean(abs(A(1:length(B))-B)/std(A(1:128))))
    % EVM_Lior_Interpolation=db(mean(abs(A(1:156)-B(1:156))/std(A(1:156))))
    %    EVM_Lior_Interpolation1=db(mean(abs(A(1:156)-D(1:156))/std(A(1:156))))
    % stem(EVM_Lior)
end

%% 12) Preamble sync enhancement
%%% I perceived that the interpolation act was by itself enhancing the
%%% preamble_synch somehow. therefore, the enhancement needs to be done
%%% after the concatenation and the interpolation, and without by mistake
%%% enhancing part of the next chips, that belong to the preamble_CE.
%%% Hence, we need to find the exact place of the samples we would like to
%%% enhance. there might be a problem since, the interpolation action "smears" the limit be
%%% tween the samples we would like to enhance and those we would not.
%%% for that purpose, i decided to make use of the 4 zeros left at
%%% the end of the premable_synch at the OFDM modulator output,
%%% Limits of the premable_synch sequence within the signal: evolution
%%% along the chain, from OFDM modulator output to Interpolator output,
%%% where we wish to enhance it and only it

% Modultor output
limit_low_low=1; % index of lowest chip
limit_high_low=N_preamble_synch*N_FFT-4+N_CP; % index of highest non-zero chip (reminder: preamble_synch initially terminates with 4 zeros)
limit_high_high=N_preamble_synch*N_FFT+N_CP+1; % index of highest zero chip still part of the preamble

% Upsampler output

limit_low_low=1; % 1st chip remains 1st
limit_high_low=(limit_high_low-1)*N_upsample+1;
limit_high_high=(limit_high_high-1)*N_upsample+1;
% LPF output

limit_low_low=limit_low_low+round(Group_delay_interp);
limit_high_low=limit_high_low+round(Group_delay_interp);
limit_high_high=limit_high_high+round(Group_delay_interp);

% chosen index: assuring that we are both above the highest limit of
% nonzero elements, and below the highest of the zero elements, after which
% the other preamble begins
limit_middle=round(mean([limit_high_low,limit_high_high]));

%%% The enhancement
P_max_preamble_synch=max(abs(Signal_Tx(limit_low_low-10:limit_middle)).^2);
P_avg_preamble_synch=mean(abs(Signal_Tx(limit_low_low-10:limit_middle)).^2);
PAPR_preamble_synch=db(P_max_preamble_synch/P_avg_preamble_synch,'power');

P_max_data=max(abs(Signal_Tx(limit_high_high:end)).^2);
P_avg_data=mean(abs(Signal_Tx(limit_high_high:end)).^2);
PAPR_data=db(P_max_data/P_avg_data,'power');


%%% power enhancement
if db(P_max_data,'power')-db(P_max_preamble_synch,'power')<Enhancement_prmbl_synch
    error('impossible to enhance preambles Synch')
end

Signal_Tx_Pre_enhancement=Signal_Tx;
Signal_Tx(limit_low_low-10:limit_middle)=Signal_Tx(limit_low_low-10:limit_middle)*db2mag(Enhancement_prmbl_synch); % the actual enhancement


if test_signal_processing_flag
    
    close all
    figure
    set(gcf,'windowstyle','docked')
    subplot(2,1,1)
    plot(db(abs(Signal_Tx_Pre_enhancement).^2,'power'));
    title('before enhancement')
    subplot(2,1,2)
    plot(db(abs(Signal_Tx).^2,'power'));
    title('after enhancement')
    
end

%%% power ratio calculations

E9=sum(sum((abs(Signal_Tx_Pre_enhancement)).^2)); %
E10=sum(sum((abs(Signal_Tx)).^2)); %

R10=E9/E10; %

%% 13) IF Up-Conversion

if strcmp(BB_IF_flag,'IF')
    
    Ts=1/Fs;
    t=0:Ts:(length(Signal_Tx)-1)*Ts;
    
    if (F_if+B_signal_net_Hz>Fs/2 || F_if-B_signal_net_Hz <0) % signal should not exceed nyquist zone to the right or dC frequencies to the left
        error('F_if and BW badly adpated to Fs after interpolation')
    end
    
    
    
    Signal_Tx=real(Signal_Tx.*exp(j*2*pi*F_if*t')); % accumulated group delay=350
    
    Signal_Tx_upconverted=Signal_Tx;
    
    %%% Start: Uncomment In case of realization over FW or SW with fixed point elements
    % Delta_F=F_chip*N_upsample/2^N_phase_NCO;% NCO's Frequency resolution. T sample at NCO is: Ts=1/(F_chip*N_upsample)
    % P=round(F_if/Delta_F);
    % F_if=P*Delta_F; % the actual F_if
    %
    % Delta_Phase=F_if*2^N_phase_NCO/(F_chip*N_upsample); % Sets the Fout see NCO's HELP
    % if (mod(Delta_Phase,1)~=0)
    %     error('NCOs phase increment is not an integer number')
    % end
    
    %%% End
    
    if test_signal_processing_flag
        % Pilots might spread along several bins becuase the bins of FFT applied over the whole
        % signal are narrower than those of the OFDM FFT, and therefor hte energy is spread along several bins
        figure(4)
        set(gcf,'windowstyle','docked');
        
        subplot(2,1,1)
        f=linspace(-Fs/2,Fs/2,length(Signal_Tx));
        A=fftshift(fft(Signal_Tx));
        plot(f/1e3,db(abs(A)))
        title(['Signal Tx after upconverting by ',num2str(F_if/1e3),' [kHz]'])
        xlabel('freq [kHz]')
        grid on
        grid minor
        
        subplot(2,1,2)
        plot(db(abs(Signal_Tx).^2,'power'));
        title('after upconversion')
        
        
    end
    
    
    % phase_offset=0;
    
    
else
    Signal_Tx_upconverted=Signal_Tx;
    
end

%% 14) Inverse Sinc Filter

if ~strcmp(Configuration,'Calibration') && ~Ref_Cfg_flag% in any case other than Calibration, the interpolation factor and thus filter need to be flexible in order to adapt to the variable F_chip and the maximum Fsample cnstrained by the NI 6212
    
    
    Fpass=F_if+F_chip/2;
    Fstop=(Fs/2)*0.99;
    
    H_inv_sinc_SPEC=fdesign.isinclp('Fp,Fst,Ap,Ast',Fpass,Fstop,0.001,55,Fs);%672
    H_inv_sinc_flt=design(H_inv_sinc_SPEC);
    
elseif strcmp(Configuration,'Calibration') % calibration
    load('Filters.mat') % load all the filters at once
    
end

%Signal_Tx=filter(H_inv_sinc_flt,Signal_Tx); % accumulated group delay=351+140=491
Signal_Tx=conv(H_inv_sinc_flt.numerator,Signal_Tx,'full'); %  % conv instead of filter so as not to loose the suffix samples of the signal, that might cause loosing an actual payload frame
Group_delay_inv_sinc=H_inv_sinc_flt.order/2;
Group_delay_total=Group_delay_total+Group_delay_inv_sinc;

%%% Power ratios
E11_1=sum(sum((abs(Signal_Tx_upconverted)).^2));
E11_2=sum(sum((abs(Signal_Tx)).^2)); % Total OFDM block (with pilots and guard bands, but w/o CP) Energy, without the preamble_CE
R11=E11_1/E11_2; %


%% 14)Differentiator Filter

if Pre_emphasis_flag% in any case other than Calibration, the interpolation factor and thus filter need to be flexible in order to adapt to the variable F_chip and the maximum Fsample cnstrained by the NI 6212

    H_diff_SPEC = fdesign.differentiator(41); % Filter order is 33.
    H_diff_flt = design(H_diff_SPEC,'firls');
    
    %fvtool(H_diff_flt,'magnitudedisplay','zero-phase')
    Signal_Tx=conv(H_diff_flt.numerator,Signal_Tx,'full'); %  % conv instead of filter so as not to loose the suffix samples of the signal, that might cause loosing an actual payload frame
    Group_delay_diff=H_diff_flt.order/2;
    Group_delay_total=Group_delay_total+Group_delay_diff;

    
   
    if test_signal_processing_flag
        
        figure(3)
        set(gcf,'windowstyle','docked');
        f=linspace(-Fs/2,Fs/2,length(Signal_Tx));
        C=fftshift(fft(Signal_Tx,length(f)));
        plot(f/1e3,db(abs(C)))
        xlabel('freq [kHz]')
        grid on
        grid minor
        
    end
    
end
%% 15) update of structures

OFDM_matrix_Tx_f_no_CP_no_Preamble=OFDM_matrix_Tx_f;


Testing_data.Group_delay_total=Group_delay_total;
Testing_data.Group_delay_interp=Group_delay_interp;
Testing_data.OFDM_matrix_Tx_f_no_CP_no_Preamble=OFDM_matrix_Tx_f_no_CP_no_Preamble;

OFDM_config.OFDM_Ref_symbol_Tx_f=OFDM_Ref_symbol_Tx_f;
OFDM_config.Pilots_matrix_Tx=Pilots_matrix_Tx;
OFDM_config.F_chip=F_chip; % calibration case F_chip might vary
OFDM_config.R1=R1;
OFDM_config.R2=R2;
OFDM_config.R3=R3;
OFDM_config.R4=R4;
OFDM_config.R5=R5;
OFDM_config.R6=R6;
OFDM_config.R7=R7;
OFDM_config.R10=R10;
OFDM_config.R11=R11;
OFDM_config.R_lengths=R_lengths;
% OFDM_config.E9=E9;

IF_chain_config.N_upsample=N_upsample;


end

