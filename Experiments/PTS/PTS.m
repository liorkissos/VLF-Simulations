
clc
clear
close all

debug=0;




%% User defined parameters

P=1000; % number of iterations

reference_cfg=1;

%scrambling='interleaved';
scrambling='contiguous';

%%% Standard. Hee& Han values: N=64, 128 . Modulation= QPSK (p.3)
N=2^6; % length of symbol
Modulation=4; % modulation index

hMod_data = comm.RectangularQAMModulator('ModulationOrder',Modulation,'NormalizationMethod','Average Power','BitInput',0);
hDemod = comm.RectangularQAMDemodulator('ModulationOrder',Modulation,'NormalizationMethod','Average Power','BitOutput',0);

%%% PTS.  Hee& Han values:  M=8, W=4
M=4; % number of blocks the symbol is divided into
W=4; % number of possible phase
L=4; % interpolation factor. see Jiang& Wu equation 5


if reference_cfg
    N=2^6; % length of symbol
    Modulation=4; % modulation index
    
    M=8; % number of blocks the symbol is divided into
    W=4; % number of possible phase
    L=4; % interpolation factor. see Jiang& Wu equation 5
    
end

%% Sanity checks

if mod(log2(M),1) || mod(log2(N),1)
    error('M and N should be a power of 2')
end

if mod(N/M,1)
    error('the symbol s subset lengths should be an integer number')
end



Signal_orig=[];
Signal_manip=[];
for pp=1:P
    
    
    %% Data Generation
    
    data_Tx=randi(Modulation,N,1)-1;
    
    X_orig = step(hMod_data, data_Tx); % QAM modulation
    
    X_orig_shift=ifftshift(X_orig);
    
    %% PAPR calculation: Pre Processing
    
    x_Tx_orig=ifft(X_orig_shift);
   % x_Tx_orig=interp(x_Tx_orig,L);
    [PAPR_orig]=PAPR_calc(x_Tx_orig,L);
    
    if debug
        %%% Frequency domain interpolation (as was demonstrated in the articles) versus time domain interpolation
        %%% there is some difference (T domain PAPR is higher)that might be due
        %%% to the interpolation filter and the windowing. the frequency domain
        %%% zero-padding+idft is an "ideal" time domain interpolation is some
        %%% way
        x_Tx_orig1=L*ifft(X_orig_shift,L*N); % frequency domain interpolation, by zero padding
        [PAPR_orig1,aaa]=PAPR_calc(x_Tx_orig1,1)
        plot(real(x_Tx_orig1));hold on;plot(real(x_Tx_orig))
        legend('F domain interp','T domain interp')
        
    end
    
    %% PTS Processing 1:Tx Symbol splitting
    
    X_subs_mat=zeros(N,M);
    
    switch scrambling
        
        case 'interleaved'
            
            for kk=1:M %%% splitting in a "polyphase style"
                ind_vec=kk:M:N;
                X_subs_mat(:,kk)=upsample(X_orig_shift(ind_vec),M,kk-1);
                
            end
            
        case 'contiguous'
            
            for kk=1:M
                X_subs_mat(:,kk)=[zeros((kk-1)*N/M,1);X_orig_shift((kk-1)*N/M+1:(kk)*N/M);zeros((M-kk)*N/M,1)];
            end
            
    end
    
    %% PTS Processing 2: Phase optimization: see PAPR Reduction of OFDM Signals Using a Reduced Complexity PTS Technique, Hee Han& Hong Lee
    
    PAPR_min_vec=zeros(M,1);
    PAPR_min_vec(1)=PAPR_orig;
    
    b=ones(M,1); % see page 2 bottom left
    for     kk=2:M
        
        x_subs_mat=ifft(X_subs_mat); % IDFT over the zero padded subsets of X
        x_Tx_mat=x_subs_mat*diag(b); % like multiplying every column by the same integer
        x_Tx=sum(x_Tx_mat,2);
        
        PAPR_ref=PAPR_calc(x_Tx,L);
        PAPR_min=PAPR_ref;
        
        ll_opt=0;
        for ll=0:W-1 % equation (4) in Han & Lee
            
            b(kk)=exp(j*2*pi*ll/W);
            x_subs_mat=ifft(X_subs_mat);
            x_Tx_mat=x_subs_mat*diag(b);
            x_Tx=sum(x_Tx_mat,2);
            
            PAPR=PAPR_calc(x_Tx,L);
            
            if PAPR<PAPR_min
                PAPR_min=PAPR;
                ll_opt=ll;
            end
            
        end
        
        b(kk)=exp(j*2*pi*ll_opt/W);
        
        
        PAPR_min_vec(kk)=PAPR_min;
        
    end
    
    %% PTS Processing 3: IDFT and summation
    
    x_subs_mat=ifft(X_subs_mat);
    x_Tx_mat=x_subs_mat*diag(b);
    x_Tx=sum(x_Tx_mat,2);
    
    %% PAPR measurement: Post Processing
    
    PAPR_manip(pp)=PAPR_calc(x_Tx,L);
    PAPR_reduction(pp)=PAPR_orig-PAPR_manip(pp);
    
    %% Signal saving
    
    Signal_orig=[Signal_orig;x_Tx_orig];
    Signal_manip=[Signal_manip;x_Tx];
    
end

%% Analysis (of the whole signal)

ccdf = comm.CCDF('AveragePowerOutputPort',true, ...
    'PeakPowerOutputPort',true,'MaximumPowerLimit',25);

[CCDFy,CCDFx,AvgPwr,PeakPwr] = ccdf([Signal_orig Signal_manip]);

figure
plot(ccdf)
set(gcf,'windowstyle','docked')
shg
legend('Original','Manipulated')


%% Rx (of a single Tx symbol)

X_Rx=fft(x_Tx);

X_Rx_mat=zeros(N,M);

switch scrambling
    
    case 'contiguous'
        
        for kk=1:M
            X_Rx_mat(:,kk)=(1./b(kk))*[zeros((kk-1)*N/M,1);X_Rx((kk-1)*N/M+1:(kk)*N/M);zeros((M-kk)*N/M,1)];
        end
        
    case 'interleaved'
        
        for kk=1:M
            ind_vec=kk:M:N;
            X_Rx_subs=X_Rx(ind_vec);
            X_Rx_mat(:,kk)=(1./b(kk))*upsample(X_Rx_subs,M,kk-1);
            
        end
        
end

X_Rx=sum(X_Rx_mat,2);

X_Rx=fftshift(X_Rx);

data_Rx=step(hDemod,X_Rx);

 Analysis

error=max(abs(data_Tx-data_Rx));
if error
    error('Tx and Rx manipulations are not the inverse of each other')
end

figure
set(gcf,'windowstyle','docked')
shg
plot(PAPR_min_vec)
grid minor
title({['PAPR optimization evolution'],['PAPR reduction=',num2str(PAPR_reduction),'[dB]']})
ylabel('PAPR [dB]')
xlabel('iteration')


if debug
    
    figure
    set(gcf,'windowstyle','docked')
    shg
    stem(data_Tx);hold on;stem(data_Rx)
    title(['Error= ',num2str(error),''])
    legend('original','reconstructed')
    grid minor
    
end


