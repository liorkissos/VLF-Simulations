
clc
clear
close all

debug=0;


%% User defined parameters

%%% Standard
N=2^8; % length of symbol
Modulation=16; % modulation index

hMod_data = comm.RectangularQAMModulator('ModulationOrder',Modulation,'NormalizationMethod','Average Power','BitInput',0);
hDemod = comm.RectangularQAMDemodulator('ModulationOrder',Modulation,'NormalizationMethod','Average Power','BitOutput',0);

%%% PTS
M=32; % number of blocks the symbol is divided into
W=8; % number of possible phase

%% Sanity checks

if mod(log2(M),1) || mod(log2(N),1)
    error('M and N should be a power of 2')
end

if mod(N/M,1)
    error('the symbol s subset lengths should be an integer number')
end

%% Data Generation

data_Tx=randi(Modulation,N,1)-1;

X_orig = step(hMod_data, data_Tx); % QAM modulation

X_orig_shift=ifftshift(X_orig);

%%
x_Tx_orig=ifft(X_orig_shift);
P_max_orig=max(x_Tx_orig.*conj(x_Tx_orig));
P_avg_orig=mean(x_Tx_orig.*conj(x_Tx_orig));
PAPR_orig=10*log10(P_max_orig/P_avg_orig);


%% Tx Symbol splitting

%%% splitting in a "polyphase style"

X_subs_mat=zeros(N,M);
for kk=1:M
    ind_vec=kk:M:N;
    X_subs_mat(:,kk)=upsample(X_orig_shift(ind_vec),M,kk-1);
    
end

%% Phase optimization: see PAPR Reduction of OFDM Signals Using a Reduced Complexity PTS Technique, Hee Han& Hong Lee

PAPR_min_vec=zeros(M,1);
PAPR_min_vec(1)=PAPR_orig;

b=ones(M,1); % see page 2 bottom left
for     kk=2:M
    
    x_subs_mat=ifft(X_subs_mat); % IDFT over the zero padded subsets of X
    x_Tx_mat=x_subs_mat*diag(b); % like multiplying every column by the same integer
    x_Tx=sum(x_Tx_mat,2);
    
    P_max=max(x_Tx.*conj(x_Tx));
    P_avg=mean(x_Tx.*conj(x_Tx));
    PAPR_ref=10*log10(P_max/P_avg);
    PAPR_min=PAPR_ref;
    
    ll_opt=0;
    for ll=0:W-1 % equation (4) in Han & Lee
        
        b(kk)=exp(j*2*pi*ll/W);
        x_subs_mat=ifft(X_subs_mat);
        x_Tx_mat=x_subs_mat*diag(b);
        x_Tx=sum(x_Tx_mat,2);
        
        P_max=max(x_Tx.*conj(x_Tx));
        P_avg=mean(x_Tx.*conj(x_Tx));
        PAPR=10*log10(P_max/P_avg);
        
        if PAPR<PAPR_min
            PAPR_min=PAPR;
            ll_opt=ll;
        end
        
    end
    
    b(kk)=exp(j*2*pi*ll_opt/W);
    
    
    PAPR_min_vec(kk)=PAPR_min;
    
end

x_subs_mat=ifft(X_subs_mat);
x_Tx_mat=x_subs_mat*diag(b);
x_Tx=sum(x_Tx_mat,2);


%% PAPR measurement

% x_Tx_orig=ifft(X_orig_shift);
% P_max_orig=max(x_Tx_orig.*conj(x_Tx_orig));
% P_avg_orig=mean(x_Tx_orig.*conj(x_Tx_orig));
% PAPR_orig=10*log10(P_max_orig/P_avg_orig)

P_max_manip=max(x_Tx.*conj(x_Tx));
P_avg_manip=mean(x_Tx.*conj(x_Tx));
PAPR_manip=10*log10(P_max_manip/P_avg_manip);
PAPR_reduction=PAPR_orig-PAPR_manip

%% Rx

X_Rx=fft(x_Tx);

X_Rx_mat=zeros(N,M);
for kk=1:M
    ind_vec=kk:M:N;
    X_Rx_subs=X_Rx(ind_vec);
    X_Rx_mat(:,kk)=(1./b(kk))*upsample(X_Rx_subs,M,kk-1);
    
end

X_Rx=sum(X_Rx_mat,2);

X_Rx=fftshift(X_Rx);

data_Rx=step(hDemod,X_Rx);


%% Analysis

error=max(abs(data_Tx-data_Rx))

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


