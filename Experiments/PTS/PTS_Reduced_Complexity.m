
clc
clear
close all

debug=0;


%% User defined parameters

P=1000; % number of symbols in the stream

reference_cfg=1;

PTS_algorithm= 'Iterative_Flipping';
%PTS_algorithm= 'Reduced_Complexity';

%scrambling='interleaved';
scrambling='contiguous';

%%% Standard. Hee& Han values: N=64, 128 . Modulation= QPSK (p.3)
N=2^6; % length of symbol
Modulation=4; % modulation index

hMod_data = comm.RectangularQAMModulator('ModulationOrder',Modulation,'NormalizationMethod','Average Power','BitInput',0);
hDemod = comm.RectangularQAMDemodulator('ModulationOrder',Modulation,'NormalizationMethod','Average Power','BitOutput',0);

%%% PTS.  Hee& Han values:  M=8, W=4
M=4; % number of blocks the symbol is divided into
W=4; % number of possible phases
L=4; % interpolation factor. see Jiang& Wu equation 5


if reference_cfg
    N=2^6; % length of symbol
    Modulation=4; % modulation index
    
    M=8; % number of blocks the symbol is divided into
    W=4; % number of possible phases
    L=4; % interpolation factor. see Jiang& Wu equation 5
    r=2; % Hamming search radius of the Reduced Complexity PTS
    
end

%% Sanity checks

if mod(log2(M),1) || mod(log2(N),1)
    error('M and N should be a power of 2')
end

if mod(N/M,1)
    error('the symbol s subset lengths should be an integer number')
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%% Algorithm execution
%%%%%%%%%%%%%%%%%%%%%%%%

Signal_orig=[];
Signal_manip=[];

PAPR_manip=1000*ones(P,1);
%  PAPR_reduction=1000*ones(P,1);

for pp=1:P % running over the stream
    
    
    %% Data Generation
    
    data_Tx=randi(Modulation,N,1)-1;
    
    X_orig = step(hMod_data, data_Tx); % QAM modulation
    
    X_orig_shift=ifftshift(X_orig);
    
    %% PAPR calculation: Pre Processing
    
    x_Tx_orig=ifft(X_orig_shift);
    % x_Tx_orig=interp(x_Tx_orig,L);
    [PAPR_orig(pp)]=PAPR_calc(x_Tx_orig,L);
    
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
    
    x_subs_mat=ifft(X_subs_mat); % IDFT over the zero padded subsets of X
    
    
    %% PTS Processing 2: Phase optimization: see PAPR Reduction of OFDM Signals Using a Reduced Complexity PTS Technique, Hee Han& Hong Lee
    %%% the technique: 1) we divide a block into M subsets. 2) we multiply
    %%% subset # 1 with 1. 2) with each of the follwing M-1 subsets (index kk) we operate as
    %%% follows: We look for an optimum complex factor to minimize the total
    %%% PAPR with the symbol so far . we begin with a b vector that is a
    %%% "ones" vector. b(1)=1 and we do not change it. with b(2) and on, we
    %%% run over the W possible phases to choose the one that optimizes the
    %%% PAPR
    
    PAPR_min_vec=zeros(M,1);
    PAPR_min_vec(1)=PAPR_orig(pp);
    
    b=ones(M,1); % see page 2 bottom left
    b_opt=ones(M,1);
    
    switch PTS_algorithm
        
        case 'Iterative_Flipping'
            
            %   x_subs_mat=ifft(X_subs_mat); % IDFT over the zero padded subsets of X
            
            for     kk=2:M % running over the M sub-blocks
                
                %   x_subs_mat=ifft(X_subs_mat); % IDFT over the zero padded subsets of X
                x_Tx_mat=x_subs_mat*diag(b); % like multiplying every column by the same integer
                x_Tx=sum(x_Tx_mat,2);
                
                PAPR_ref=PAPR_calc(x_Tx,L);
                PAPR_min=PAPR_ref;
                
                ll_opt=0;
                for ll=0:W-1 % equation (4) in Han & Lee. running over the W possible phases
                    
                    b(kk)=exp(j*2*pi*ll/W);
                    % x_subs_mat=ifft(X_subs_mat);
                    x_Tx_mat=x_subs_mat*diag(b); % multiplication of every column of x_Tx_mat by a scalar; a term in  vector
                    x_Tx=sum(x_Tx_mat,2);
                    
                    PAPR=PAPR_calc(x_Tx,L);
                    
                    if PAPR<PAPR_min
                        PAPR_min=PAPR;
                        ll_opt=ll;
                        b_opt(kk)=exp(j*2*pi*ll_opt/W);
                    end
                    
                end
                
                b(kk)=exp(j*2*pi*ll_opt/W);
                
                
                PAPR_min_vec(kk)=PAPR_min;
                
            end
            
        case 'Reduced_Complexity'
            
            PAPR_min=PAPR_orig(pp); % PAPR is referenced to the one obtained with no signal processing
            
            ind_mat=combnk(1:M,r);
            
            for ii=1:size(ind_mat,1) % running on the matrix rows
                
                a=exp(j*2*pi*(0:W-1)/W);
                c=combnk(a,r);
                
                for jj=1:length(c)
                    
                    b(ind_mat(ii,:))=c(jj,:);
                    
                    x_Tx_mat=x_subs_mat*diag(b); % multiplication of every column of x_Tx_mat by a scalar; a term in  vector
                    x_Tx=sum(x_Tx_mat,2);
                    
                    PAPR=PAPR_calc(x_Tx,L);
                    
                    if PAPR<PAPR_min
                        PAPR_min=PAPR
                     %   ind_opt=[ii;jj];      
                        b_opt=b;
                    end
                    
                end % for: b indexes combinations
            end % for: b values

    end % Switch PTS algorithms
    
    %% PTS Processing 3: IDFT and summation
    
    x_subs_mat=ifft(X_subs_mat);
    x_Tx_mat=x_subs_mat*diag(b_opt);
    x_Tx=sum(x_Tx_mat,2);
    
    %% PAPR measurement: Post Processing
    
    PAPR_manip(pp)=PAPR_calc(x_Tx,L);
    PAPR_reduction(pp)=PAPR_orig(pp)-PAPR_manip(pp);
    
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
legend('Original','Manipulated')

figure
set(gcf,'windowstyle','docked')
plot(PAPR_orig)
hold on
plot(PAPR_manip)
legend('Pre-PTS','Post-PTS')
grid minor
title({['PAPR post PTS of a stream of N=',num2str(N),' symbols'],['average PAPR reduction=',num2str(db(mean(db2mag(PAPR_reduction)))),'[dB]']})
ylabel('PAPR [dB]')
xlabel('iteration')



%% Rx (of a single Tx symbol)

X_Rx=fft(x_Tx);

X_Rx_mat=zeros(N,M);

switch scrambling
    
    case 'contiguous'
        
        for kk=1:M
            X_Rx_mat(:,kk)=(1./b_opt(kk))*[zeros((kk-1)*N/M,1);X_Rx((kk-1)*N/M+1:(kk)*N/M);zeros((M-kk)*N/M,1)];
        end
        
    case 'interleaved'
        
        for kk=1:M
            ind_vec=kk:M:N;
            X_Rx_subs=X_Rx(ind_vec);
            X_Rx_mat(:,kk)=(1./b_opt(kk))*upsample(X_Rx_subs,M,kk-1);
            
        end
        
end

X_Rx=sum(X_Rx_mat,2);

X_Rx=fftshift(X_Rx);

data_Rx=step(hDemod,X_Rx);

%% Analysis (of a single symbol, the last one in the P long stream)

hEVM=comm.EVM('Normalization','Average reference signal power');
EVM=step(hEVM,X_orig,X_Rx);
EVM_dB=db(EVM/100);

figure
scatterplot(X_Rx);
set(gcf,'windowstyle','docked');
title(gca,['scatter plot of the last symbol of the stream. EVM=',num2str(EVM_dB),'[dB]'])
grid on
grid minor

error=max(abs(data_Tx-data_Rx));
if error
    error('Tx and Rx manipulations are not the inverse of each other')
end

figure
set(gcf,'windowstyle','docked')
plot(PAPR_min_vec)
grid minor
title(['PAPR optimization evolution (of the last symbol)'])
%title({['PAPR optimization evolution (of the last symbol)'],['PAPR reduction=',num2str(PAPR_reduction),'[dB]']})
ylabel('PAPR [dB]')
xlabel('iteration')
shg


if debug
    
    figure
    set(gcf,'windowstyle','docked')
    shg
    stem(data_Tx);hold on;stem(data_Rx)
    title(['Error= ',num2str(error),''])
    legend('original','reconstructed')
    grid minor
    
end


