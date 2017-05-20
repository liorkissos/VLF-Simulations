clear
clc
close all


% N=(2^m)-1 length PN series, constructed by a prime
% polynomial (CCF architecture [canonic representation]), and thereofre is
% a m-length sequence with exceptional Autotcorrelation properties
% see: www.gaussianwaves.com/2010/09/maximum-length-sequences-m-sequences-2
% good correlation properties are needed: 1) speed of convergence of
% algorithm (BFB p. 147, p.155. BLM p. 437) 2) steady state (execess error) value (BFB., p.153)
% it conserves the properties of orthogonality under BPSK modulation

N_FFT=512;
N_GB=2*28;
N_preamble_symbols=N_FFT-2*N_GB;

%%% generating the basic subframe

m=log2(N_FFT); % Minn & Zeng. L=N_FFT/4


%N_preamble_synch=8;
%m=log2(N_preamble_synch*N_FFT/4); % Minn & Zeng. L=N_FFT/4

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
        
    otherwise
        error('no polynomial matches required premabole length')
end

%N_preamble_symbols=(2^m-1)/log2(2); % preamble is always BPSK, thus log2(2)
hMod_preamble=comm.BPSKModulator;

PAPR_preamble_CE_min=100;
n=1;
tic
while(PAPR_preamble_CE_min>4 && toc<30)
    
    initial_conditions = randi([0 1],m,1);
    hpn = comm.PNSequence('Polynomial',Polynomial,'SamplesPerFrame', N_preamble_symbols,'InitialConditionsSource','Input port');
    
    preamble_Tx = step(hpn,initial_conditions);
    
    Preamble_stream_Tx=step(hMod_preamble,preamble_Tx);
    Preamble_stream_Tx=Preamble_stream_Tx';
    
    %%% PAPR Calc
    index_middle=(N_FFT-2*N_GB)/2+1;
    Preamble_stream_Tx_f=[zeros(1,N_GB),Preamble_stream_Tx(1:index_middle-1),0,Preamble_stream_Tx(index_middle:end),zeros(1,N_GB-1)];
    
    Preamble_stream_Tx_f=ifftshift(Preamble_stream_Tx_f,1); %  fftshift ove the columns. see comment above
    Preamble_stream_Tx_t=ifft(Preamble_stream_Tx_f); % the famous OFD IDFT operation
    
    P_max_preamble_CE=max(conj(Preamble_stream_Tx_t).*Preamble_stream_Tx_t);
    P_avg_preamble_CE=mean(conj(Preamble_stream_Tx_t).*Preamble_stream_Tx_t);
    PAPR_preamble_CE=db(P_max_preamble_CE/P_avg_preamble_CE,'power');
    
    if PAPR_preamble_CE<PAPR_preamble_CE_min
        PAPR_preamble_CE_min=PAPR_preamble_CE;
        AAA=[Preamble_stream_Tx(1:index_middle-1),0,Preamble_stream_Tx(index_middle:end)];
    end
    
    PAPR_preamble_CE_min=min(PAPR_preamble_CE_min,PAPR_preamble_CE);
    n=n+1;
end
toc


save PN_seq_512.mat AAA

%% Display
figure
set(gcf,'windowstyle','docked')
plot(xcorr(real(Preamble_stream_Tx)))
%plot(xcorr(real(AAA)))
