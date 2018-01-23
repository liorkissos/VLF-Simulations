
clc
clear
close all

debug=0;


%% User defined parameters

P=100; % number of OFDM symbols ( super symbols) in the stream

reference_cfg=1;

%PTS_Algorithm= 'Iterative Flipping';
PTS_Algorithm= 'Reduced Complexity-mine';
%PTS_Algorithm= 'Reduced Complexity-Article';

%scrambling='interleaved';
scrambling='contiguous';

%%% Standard. Hee& Han values: N=64, 128 . Modulation= QPSK (p.3)
N=2^6; % length of symbol
Modulation=4; % modulation index

hMod_data = comm.RectangularQAMModulator('ModulationOrder',Modulation,'NormalizationMethod','Average Power','BitInput',0);
hDemod = comm.RectangularQAMDemodulator('ModulationOrder',Modulation,'NormalizationMethod','Average Power','BitOutput',0);

%%% PTS.  Hee& Han values:  M=8, W=4
M=8; % number of blocks the symbol is divided into
W=4; % number of possible phases
L=4; % interpolation factor. see Jiang& Wu equation 5
r=2;

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


Signal_orig=zeros(N*P,1);
Signal_manip=zeros(N*P,1);

PAPR_orig=zeros(P,1);

PAPR_manip=1000*ones(P,1);
%  PAPR_reduction=1000*ones(P,1);

PAPR_min_stat_mat=zeros(1e6,1);

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
    
    x_subs_mat=ifft(X_subs_mat); % IDFT over the zero padded subsets of X. see Han& Lee fig. 1
    
    
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
    
    %PAPR_min_stat_vec(1)
    
    b=ones(M,1); % The phase factor vector. see page 2 bottom left.
    b_opt=ones(M,1);
    
    switch PTS_Algorithm
        
        case 'Iterative Flipping'
            
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
            
        case 'Reduced Complexity-mine' 
            
            kk=1;
            
            PAPR_min=PAPR_orig(pp); % PAPR is referenced to the one obtained with no signal processing
            
            ind_mat=combnk(1:M,r);% we check the PAPR when r terms within b vector change. ind_mat contains the possible indexes of those terms
            
            
            b=ones(M,1); % at each symbol we begin with the ones vector as PTS coefficients reference vector
            b_opt=ones(M,1);
            
            
            c=combnk(repmat(exp(j*2*pi*(0:W-1)/W),1,2),r); % each of the r terms within b has one of the values within exp(j*2*pi*(0:W-1)/W), we repmat them to enable combination of 2 identical terms, such as 1,1
            c=sort(c,2); %sorting is a necessary condition for performing unique function
            c=unique(c,'rows'); % we remove the duplicate terms in c
            
            %%% The following stacked for loops are in fact
            
            
            for ii=1:size(ind_mat,1) % running along  the matrix columns; the different indexes of b terms
                
            
                
                for jj=1:size(c,1) % running across all the b vector possibilities that are r Hamming distant from initial b (b=ones(M,1)). r Hamming distant between 2 vectors= r terms different between the 2 vectors
                    
                    
                    if  isequal(b(ind_mat(ii,:)),c(jj,:).')  % we skip the cases where new phase factors are identical to the current ones
                        % c(jj,:)
                        continue
                    end
                    
                    
                    b(ind_mat(ii,:))=c(jj,:); % the r terms within b are assigned r values. we thus create a new vector b, different from original in r terms; a Hamming distance of r
                    
                    x_Tx_mat=x_subs_mat*diag(b); % multiplication of every column of x_Tx_mat by a scalar; a term in  vector
                    x_Tx=sum(x_Tx_mat,2);
                    
                    PAPR=PAPR_calc(x_Tx,L);
                    
                    if PAPR<PAPR_min % update the PAPR_min to the new one and move with the gradient towards the descending direction
                        PAPR_min=PAPR;
                        b_opt=b;
                    else
                        b=b_opt; % the difference from the article method: if the change is not in the gradient descet direction, stay where the last known optimum is
                    end
                    
                    
                    PAPR_min_stat_vec((ii-1)*size(c,1)+kk)=PAPR_min;
                    kk=kk+1;
                    
                end % for: b indexes combinations
                
                    kk=0;
            end % for: b values
            
            
            
            %%%% second iteration (reminder r2/I2 option, namely I2= 2 iterations). initial b is the one obtained on
            %%%% previous iteration
            
            JJ=length(PAPR_min_stat_vec);
            PAPR_min_stat_vec=PAPR_min_stat_vec';
            
            
            b=b_opt;  % initial b vecotr on 2nd iteration is the optimum of the previous one
            
            for ii=1:size(ind_mat,1) % running on the matrix rows
                
               
                
                for jj=1:size(c,1)
                    
                    
                    if  isequal(b(ind_mat(ii,:)),c(jj,:).')  % we skip the cases where new phase factors are identical to the current ones
                        continue
                    end
                    
                    b(ind_mat(ii,:))=c(jj,:);
                    
                    x_Tx_mat=x_subs_mat*diag(b); % multiplication of every column of x_Tx_mat by a scalar; a term in  vector
                    x_Tx=sum(x_Tx_mat,2);
                    
                    PAPR=PAPR_calc(x_Tx,L);
                    
                    if PAPR<PAPR_min
                        PAPR_min=PAPR;
                        b_opt=b;
                    else
                        b=b_opt; %  the difference from the article method: if the change is not in the gradient descet direction, stay where the last known optimum is
                    end
                    
                    PAPR_min_stat_vec(JJ+(ii-1)*size(c,1)+kk)=PAPR_min;
                    kk=kk+1;
                %    PAPR_min_stat_vec(JJ+ii+jj-1)=PAPR_min;
                    
                end % for: b indexes combinations
                
                 kk=0;
            end % for: b values
            
            PAPR_min_stat_vec(PAPR_min_stat_vec==0)=[];
            PAPR_min_stat_mat=PAPR_min_stat_mat(1:min(length(PAPR_min_stat_vec),size(PAPR_min_stat_mat,1)),:);
            PAPR_min_stat_vec=PAPR_min_stat_vec(1:size(PAPR_min_stat_mat,1));
            PAPR_min_stat_mat(:,pp)=PAPR_min_stat_vec;
            PAPR_min_stat_vec=[];
            
        case 'Reduced Complexity-Article' % Han& Lee, section C, page 3 upper left. r2/I2 option
            
            PAPR_min=PAPR_orig(pp); % PAPR is referenced to the one obtained with no signal processing
            
            ind_mat=combnk(1:M,r);% we check the PAPR when r terms within b vector change. ind_mat contains the possible indexes of those terms
            
            
            b=ones(M,1); % at each symbol we begin with the ones vector as PTS coefficients reference vector
            b_opt=ones(M,1);
            
            
            c=combnk(repmat(exp(j*2*pi*(0:W-1)/W),1,2),r); % each of the r terms within b has one of the values within exp(j*2*pi*(0:W-1)/W), we repmat them to enable combination of 2 identical terms, such as 1,1
            c=sort(c,2); %sorting is a necessary condition for performing unique function
            c=unique(c,'rows'); % we remove the duplicate terms in c
            
            %%% The following stacked for loops are in fact
            
            
            for ii=1:size(ind_mat,1) % running along  the matrix columns; the different indexes of b terms
                
                for jj=1:size(c,1) % running across all the b vector possibilities that are r Hamming distant from initial b (b=ones(M,1)). r Hamming distant between 2 vectors= r terms different between the 2 vectors
                    
                    b=ones(M,1); % we return to 0 at every iteration
                    
                    if  isequal(b(ind_mat(ii,:)),c(jj,:).')  % we skip the cases where new phase factors are identical to the current ones
                        % c(jj,:)
                        continue
                    end

                    b(ind_mat(ii,:))=c(jj,:); % the r terms within b are assigned r values. we thus create a new vector b, different from original in r terms; a Hamming distance of r
                    
                    x_Tx_mat=x_subs_mat*diag(b); % multiplication of every column of x_Tx_mat by a scalar; a term in  vector
                    x_Tx=sum(x_Tx_mat,2);
                    
                    PAPR=PAPR_calc(x_Tx,L);
                    
                    if PAPR<PAPR_min % update the PAPR_min to the new one and move with the gradient towards the descending direction
                        PAPR_min=PAPR;
                        b_opt=b;
%                     else
%                         b=b_opt; % if the change is not in the gradient descet direction, stay where the last known optimum is
                    end
                    
                end % for: b indexes combinations
            end % for: b values
            
            
            
            %%%% second iteration (reminder r2/I2 option, namely I2= 2 iterations). initial b is the one obtained on
            %%%% previous iteration
            
            b_opt_1st=b_opt;
            
            b=b_opt_1st;  % initial b vecotr on 2nd iteration is the optimum of the previous one
            
            for ii=1:size(ind_mat,1) % running on the matrix rows
               
                for jj=1:size(c,1)
                    
                    b=b_opt_1st;
                    
                    if  isequal(b(ind_mat(ii,:)),c(jj,:).')  % we skip the cases where new phase factors are identical to the current ones
                        continue
                    end
                    
                    b(ind_mat(ii,:))=c(jj,:);
                    
                    x_Tx_mat=x_subs_mat*diag(b); % multiplication of every column of x_Tx_mat by a scalar; a term in  vector
                    x_Tx=sum(x_Tx_mat,2);
                    
                    PAPR=PAPR_calc(x_Tx,L);
                    
                    if PAPR<PAPR_min
                        PAPR_min=PAPR;
                        b_opt=b;
%                     else
%                         b=b_opt;
                    end
                    
                end % for: b indexes combinations
            end % for: b values
            
    end % Switch PTS algorithms
    
    %% PTS Processing 3: multiplication by PTS factors and summation
    
    
    x_Tx_mat=x_subs_mat*diag(b_opt);
    x_Tx=sum(x_Tx_mat,2);
    
    %% PAPR measurement: Post Processing
    
    PAPR_manip(pp)=PAPR_calc(x_Tx,4);
    PAPR_reduction(pp)=PAPR_orig(pp)-PAPR_manip(pp);
    
    %% Signal saving
    
    Signal_orig((pp-1)*N+1:pp*N)=x_Tx_orig;
    Signal_manip((pp-1)*N+1:pp*N)=x_Tx;
       
end

%% Analysis (of the whole signal)

ccdf = comm.CCDF('AveragePowerOutputPort',true, ...
    'PeakPowerOutputPort',true,'MaximumPowerLimit',25);

[CCDFy,CCDFx,AvgPwr,PeakPwr] = ccdf([Signal_orig Signal_manip]);



figure
plot(ccdf)
set(gcf,'windowstyle','docked')
legend('Original','Manipulated')
if  strcmp(PTS_Algorithm,'Iterative Flipping')
    title({['CCDF curve. Algorithm= ',PTS_Algorithm,' '],...
        ['Modulation=QPSK'],....
        ['N fft=',num2str(N),''],...
        })
else
    title({['CCDF curve. Algorithm= ',PTS_Algorithm,' '],...
        ['Modulation=QPSK '],...
        ['N fft=',num2str(N),'. M=',num2str(M),'. W=',num2str(W),'. r=',num2str(r),'. I=2']})
end

xlabel('PAPR0 [dB]')
ylabel('Pr(PAPR>PAPR0)')



figure
set(gcf,'windowstyle','docked')
plot(PAPR_orig)
hold on
plot(PAPR_manip)
legend('Pre-PTS','Post-PTS')
grid minor
if  strcmp(PTS_Algorithm,'Iterative Flipping')
    title({['PAPR Vs. OFDM Symbol Index. Algorithm= ',PTS_Algorithm,' '],...
        ['Modulation=QPSK'],....
        ['N fft=',num2str(N),'']})
    %    ['average PAPR reduction=',num2str(db(mean(db2pow(PAPR_reduction)),'power')),'[dB]'] })
else
    title({['PAPR Vs. OFDM Symbol Index. Algorithm= ',PTS_Algorithm,' '],...
        ['Modulation=QPSK'],...
        ['N fft=',num2str(N),'. M=',num2str(M),'. W=',num2str(W),'. r=',num2str(r),'. I=2']})
    %    ['average PAPR reduction=',num2str(db(mean(db2pow(PAPR_reduction)),'power')),'[dB]']})
end
%title({['PAPR post PTS of a stream- of N=',num2str(N),' symbols'],...
%    ['average PAPR reduction=',num2str(db(mean(db2pow(PAPR_reduction)),'power')),'[dB]']})
ylabel('PAPR [dB]')
xlabel('Symbol #')

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


%figure
scatterplot(X_Rx);
set(gcf,'windowstyle','docked');
title(gca,['scatter plot of the last symbol of the stream. EVM=',num2str(EVM_dB),'[dB]'])
grid on
grid minor


error=max(abs(data_Tx-data_Rx));
if error
    error('Tx and Rx manipulations are not the inverse of each other')
end

if strcmp(PTS_Algorithm,'Iterative Flipping')
    
    figure
    set(gcf,'windowstyle','docked')
    plot(PAPR_min_vec)
    grid minor
    title(['PAPR optimization evolution (of the last symbol)'])
    %title({['PAPR optimization evolution (of the last symbol)'],['PAPR reduction=',num2str(PAPR_reduction),'[dB]']})
    ylabel('PAPR [dB]')
    xlabel('iteration')
    
end

if strcmp(PTS_Algorithm, 'Reduced Complexity-mine')
    plot(PAPR_min_stat_mat)
    set(gcf,'windowstyle','docked')
    grid minor
    title('evolution of PAPR with algorithm iterations')
    ylabel('PAR [dB]')
    xlabel('Iteration number')
end

if debug
    
    figure
    set(gcf,'windowstyle','docked')
    shg
    stem(data_Tx);hold on;stem(data_Rx)
    title(['Error= ',num2str(error),''])
    legend('original','reconstructed')
    grid minor
    
end


