function [ OFDM_matrix_PTS_Tx_t, b_opt_mat ] = PTS_Tx( OFDM_matrix_Tx_f,OFDM_config )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%% PTS Processing 2: Phase optimization: see PAPR Reduction of OFDM Signals Using a Reduced Complexity PTS Technique, Hee Han& Hong Lee
%%% the technique: 1) we divide a block into M subsets. 2) we multiply
%%% subset # 1 with 1. 2) with each of the follwing M-1 subsets (index kk) we operate as
%%% follows: We look for an optimum complex factor to minimize the total
%%% PAPR with the symbol so far . we begin with a b vector that is a
%%% "ones" vector. b(1)=1 and we do not change it. with b(2) and on, we
%%% run over the W possible phases to choose the one that optimizes the
%%% PAPR


debug=0;

%P=1000; % number of symbols in the stream
P=size(OFDM_matrix_Tx_f,2); % number of OFDM symbols in the stream=columns in the OFDM matrix
N=size(OFDM_matrix_Tx_f,1); % number of QAM symbols in the OFDM symbol=rows in the OFDM matrix

r=2; % Hamming search radius of the Reduced Complexity PTS

reference_cfg=0;



%PTS_algorithm= 'Iterative_Flipping';
%PTS_algorithm= 'Reduced_Complexity';

%scrambling='interleaved';
%scrambling='contiguous';

%%% Standard. Hee& Han values: N=64, 128 . Modulation= QPSK (p.3)
%N=2^6; % length of symbol
%Modulation=4; % modulation index

% hMod_data = comm.RectangularQAMModulator('ModulationOrder',Modulation,'NormalizationMethod','Average Power','BitInput',0);
% hDemod = comm.RectangularQAMDemodulator('ModulationOrder',Modulation,'NormalizationMethod','Average Power','BitOutput',0);

%%% PTS.  Hee& Han values:  M=8, W=4
% M=4; % number of blocks the symbol is divided into
% W=4; % number of possible phases
% L=4; % interpolation factor. see Jiang& Wu equation 5

%% User defined params

PTS_Algorithm=OFDM_config.PTS.PTS_Algorithm;
scrambling=OFDM_config.PTS.scrambling;

M=OFDM_config.PTS.M_PTS; % number of blocks the symbol is divided into
W=OFDM_config.PTS.W_PTS; % number of possible phases
L=OFDM_config.PTS.L_PTS; % interpolation factor. see Jiang& Wu equation 5

P_data=OFDM_config.PTS.P_data; %required power from
N_preamble_CE=OFDM_config.N_preamble_CE;

if reference_cfg
    N=2^6; % length of symbol
    Modulation=4; % modulation index
    
    M=8; % number of blocks the symbol is divided into
    W=4; % number of possible phases
    L=4; % interpolation factor. see Jiang& Wu equation 5
    r=2; % Hamming search radius of the Reduced Complexity PTS
    
end

if debug & 0
    
    N=256;
    M_modulation=16;
    
    M=16;
    
    hMod_data = comm.RectangularQAMModulator('ModulationOrder',M_modulation,'NormalizationMethod','Average Power','BitInput',0);
    
    OFDM_matrix_Tx_f=zeros(N,P);
end

%% Sanity checks

if mod(log2(M),1) || mod(log2(N),1)
    warning('M and N should be a power of 2')
end

if mod(N/M,1)
    error('the symbol s subset lengths should be an integer number')
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%% Algorithm execution
%%%%%%%%%%%%%%%%%%%%%%%%

%% Phase 1: PTS coeffiecients calculation

Signal_orig=[];
Signal_manip=[];

%%% initializing the output matrix: it is a time domain matrix,in which we do not want to modify the N_preamble_CE preambles.
%%% thus, we apply ifft to these symbols, and concatenate them to a zero
%%% matrix.
OFDM_matrix_PTS_Tx_t=cat(2,ifft(OFDM_matrix_Tx_f(:,1:N_preamble_CE),[],1),zeros(N,P-2)); %

b_opt_mat=zeros(M,P-N_preamble_CE);

PAPR_manip=1000*ones(P,1);

PAPR_orig=zeros(P,1);

for pp=N_preamble_CE+1:P % running over the Tx matrix, excluding the 1st& 2nd the Preamble CE
    
    
    %% Data Generation
    
    X_orig=OFDM_matrix_Tx_f(:,pp); % the entire matrix is ifftshifted before entering the function
    
    
    if debug & 0 % enable only when wanting to generate a synthetic signal. e.g; when wanting to test a different modulation or a different symbol length
        
        data_Tx=randi(M_modulation,N,1)-1;
        
        X_orig = step(hMod_data, data_Tx); % QAM modulation
        
        OFDM_matrix_Tx_f(:,pp)=X_orig;
        
        X_orig=ifftshift(X_orig);
        
        %     else
        %
        %         X_orig=OFDM_matrix_Tx_f(:,pp); % the entire matrix is ifftshifted before entering the function
        
    end
    %% PAPR calculation: Pre Processing
    
    x_Tx_orig=ifft(X_orig);
    
    [PAPR_orig(pp)]=PAPR_calc(x_Tx_orig,L);
    
    if debug & 0
        %%% Frequency domain interpolation (as was demonstrated in the articles) versus time domain interpolation
        %%% there is some difference (T domain PAPR is higher)that might be due
        %%% to the interpolation filter and the windowing. the frequency domain
        %%% zero-padding+idft is an "ideal" time domain interpolation is some
        %%% way
        x_Tx_orig1=L*ifft(X_orig,L*N); % frequency domain interpolation, by zero padding
        [PAPR_orig1,aaa]=PAPR_calc(x_Tx_orig1,1);
        plot(real(x_Tx_orig1));hold on;plot(real(x_Tx_orig))
        legend('F domain interp','T domain interp')
        
    end
    
    %% PTS Processing 1:Tx Symbol division into sub blocks
    
    X_subs_mat=zeros(N,M);
    
    switch scrambling
        
        case 'interleaved'
            
            for kk=1:M %%% splitting in a "polyphase style"
                ind_vec=kk:M:N;
                X_subs_mat(:,kk)=upsample(X_orig(ind_vec),M,kk-1);
                
            end
            
        case 'contiguous'
            
            for kk=1:M
                X_subs_mat(:,kk)=[zeros((kk-1)*N/M,1);X_orig((kk-1)*N/M+1:(kk)*N/M);zeros((M-kk)*N/M,1)];
            end
            
    end
    
    
    %% PTS Processing 2: IDFT upon the sub-blocks
    
    x_subs_mat=ifft(X_subs_mat); % IDFT over the zero padded subsets of X. see Han& Lee fig. 1
    
    
    %% PTS Processing 3: Phase optimization: see PAPR Reduction of OFDM Signals Using a Reduced Complexity PTS Technique, Hee Han& Hong Lee
    %%% the technique: 1) we divide a block into M subsets. 2) we multiply
    %%% subset # 1 with 1. 2) with each of the follwing M-1 subsets (index kk) we operate as
    %%% follows: We look for an optimum complex factor to minimize the total
    %%% PAPR with the symbol so far . we begin with a b vector that is a
    %%% "ones" vector. b(1)=1 and we do not change it. with b(2) and on, we
    %%% run over the W possible phases to choose the one that optimizes the
    %%% PAPR
    
    PAPR_min_vec=zeros(M,1);
    PAPR_min_vec(1)=PAPR_orig(pp);
    
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
            
        case 'Reduced Complexity-mine' % Han& Lee, section C, page 3 upper left. r2/I2 option
            
            PAPR_min=PAPR_orig(pp); % PAPR is referenced to the one obtained with no signal processing
            
            ind_mat=combnk(1:M,r);% we check the PAPR when r terms within b vector change. ind_mat contains the possible indexes of those terms
            
            b=ones(M,1); % at each symbol we begin with the ones vector as PTS coefficients reference vector
            b_opt=ones(M,1);
            
            
            c=combnk(repmat(exp(j*2*pi*(0:W-1)/W),1,2),r); % each of the r terms within b has one of the values within exp(j*2*pi*(0:W-1)/W), we repmat them to enable combination of 2 identical terms, such as 1,1
            c=sort(c,2); %sorting is a necessary condition for performing unique function
            c=unique(c,'rows'); % we remove the duplicate terms in c
            
            for ii=1:size(ind_mat,1) % running on the matrix rows
                
                %    b=ones(M,1);
                %   b=b_opt;
                %  c=combnk(exp(j*2*pi*(0:W-1)/W),r); % each of the r terms within b has one of the values within exp(j*2*pi*(0:W-1)/W)
                
                for jj=1:length(c)
                    
                    if  isequal(b(ind_mat(ii,:)),c(jj,:).')  % we skip the cases where new phase factors are identical to the current ones
                        % c(jj,:)
                        continue
                    end
                    
                    b(ind_mat(ii,:))=c(jj,:); % the r terms within b are assigned r values
                    
                    x_Tx_mat=x_subs_mat*diag(b); % multiplication of every column of x_Tx_mat by a scalar; a term in  vector
                    x_Tx=sum(x_Tx_mat,2);
                    
                    PAPR=PAPR_calc(x_Tx,L);
                    
                    if PAPR<PAPR_min % update the PAPR_min to the new one and move with the gradient towards the descending direction
                        PAPR_min=PAPR;
                        b_opt=b;
                    else
                        b=b_opt; % if the change is not in the gradient descet direction, stay where the last known optimum is
                    end
                    
                end % for: b indexes combinations
            end % for: b values
            
            
            
            %%%% second iteration (reminder r2/I2 option, namely I2= 2 iterations). initial b is the one obtained on
            %%%% previous iteration
            
            %b_1st_iter=b_opt; % initial b vecotr on 2nd iteration is the optimum of the previous one
            
            b=b_opt;  % initial b vecotr on 2nd iteration is the optimum of the previous one
            
            for ii=1:size(ind_mat,1) % running on the matrix rows
                
                %                 b=b_1st_iter;
                %                 c=combnk(exp(j*2*pi*(0:W-1)/W),r);
                
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
                        b=b_opt;
                    end
                    
                end % for: b indexes combinations
            end % for: b values
            
    end % Switch PTS algorithms
    
    
    %% PTS Processing 4: Coeffieicients multiplication placing the PTS coefficients in their place
    %     %%%% the chosen places are the ones are the exteremeties of the symbol (before ifftshift)
    %     %%%% after ifftshift, they are placed elsewhere (the right one). those
    %     %%%% M places normally belong to data subcarriers, so in fact we
    %     %%%% overwrite M data subcarriers
    %
    %     %%% Do not Erase!! needed in case we want to
    %     b_opt=b_opt*(sqrt(P_data)/sqrt(mean(b_opt.*conj(b_opt)))); %  enforcement of the required average power. derived from the P_total parameter
    %
    %         X_intermediate=X_orig;
    %
    %         ind_left=OFDM_config.N_FFT/2-(OFDM_config.Nguard_band_right+1); % index of the 1st data subcarrier left of the pilot after ifftshift
    %         ind_right=ind_left+1+OFDM_config.Nguard_band_left+OFDM_config.Nguard_band_right+1+1; % index of the 1st data subcarrier right of the pilot after ifftshift
    %
    %         X_intermediate(ind_left-M/2+1:ind_left)=b_opt(1:M/2);
    %         X_intermediate(ind_right:ind_right+M/2-1)=b_opt(M/2+1:M);
    %
    %         OFDM_matrix_intermediate_Tx_f(:,pp)=X_intermediate;
    %
    %      %   error('does not work well. b_opt coeffs need to be added differently. maybe another symbol')
    
    
    %% PTS Processing 5: Coeffieicients multiplication
    
    x_Tx_mat=x_subs_mat*diag(b_opt);
    OFDM_matrix_PTS_Tx_t(:,pp)=sum(x_Tx_mat,2);
    
    b_opt_mat(:,pp-N_preamble_CE)=b_opt; %saving the coefficients for  their applicaion later on. Doing that before scaling the power so as to not affect the power of the signal, since b_opt terms have a power of 1
    
    
    if debug
        
        %%% PAPR measurement: Post Processing
        
        PAPR_manip(pp)=PAPR_calc(x_Tx,L);
        PAPR_reduction(pp)=PAPR_orig(pp)-PAPR_manip(pp);
        
        %%% Signal saving
        
        Signal_orig=[Signal_orig;x_Tx_orig];
        Signal_manip=[Signal_manip;x_Tx];
        
        
        
    end
    
    PAPR_min_stat_vec(pp-N_preamble_CE)=PAPR_min;
    
end % for running over symbols


% %%%% TEMP
% % OFDM_matrix_intermediate_Tx_f=OFDM_matrix_Tx_f;
% % error('temporary!!!!, line above needs to be in comment ')
% %%%%%%%
%
% %% Phase 2: PTS application
%
% for pp=2:P % running over the Tx matrix,excluding the Preamble CE
%
%
%     %% Data Generation
%
%     X_PTS=OFDM_matrix_intermediate_Tx_f(:,pp); % the entire matrix is ifftshifted before entering the function
%
%
%
%     %% PTS Processing :Tx Symbol splitting and multiplying by PTS coefficient
%
%     X_PTS_subs_mat=zeros(N,M);
%
%     switch scrambling
%
%         case 'interleaved'
%
%             for kk=1:M %%% splitting in a "polyphase style"
%                 ind_vec=kk:M:N;
%                 X_PTS_subs_mat(:,kk)=upsample(X_PTS(ind_vec),M,kk-1);
%                 X_PTS_subs_mat(:,kk)=X_PTS_subs_mat(:,kk)*b_opt_mat(kk,pp); % multiplying by the PTS coefficients
%             end
%
%         case 'contiguous'
%
%             for kk=1:M
%                 X_PTS_subs_mat(:,kk)=[zeros((kk-1)*N/M,1);X_PTS((kk-1)*N/M+1:(kk)*N/M);zeros((M-kk)*N/M,1)];
%                 X_PTS_subs_mat(:,kk)=X_PTS_subs_mat(:,kk)*b_opt_mat(kk,pp); % multiplying by the PTS coefficients
%             end
%
%     end
%
%     x_PTS_subs_mat=ifft(X_PTS_subs_mat); % IDFT over the zero padded subsets of X. see Han& Lee fig. 1
%
%     OFDM_matrix_PTS_Tx_t(:,pp)=sum(x_PTS_subs_mat,2); %summation along the row dimension. see Han& Lee
%
%
% end

%% Debug

% PAPR_min_stat_vec_new=PAPR_min_stat_vec;
% histogram(PAPR_min_stat_vec,100)
% load('PAPR_min_stat_vec')
%
% ecdf(PAPR_min_stat_vec)
% hold on
% ecdf(PAPR_min_stat_vec_new)

if debug
    
    Signal_manip=OFDM_matrix_PTS_Tx_t(:);
    
    A=ifft(OFDM_matrix_Tx_f);
    Signal_orig=A(:);
    
    ccdf = comm.CCDF('AveragePowerOutputPort',true, ...
        'PeakPowerOutputPort',true,'MaximumPowerLimit',10);
    
    [CCDFy,CCDFx,AvgPwr,PeakPwr] = ccdf([Signal_orig Signal_manip]);
    
    figure
    plot(ccdf)
    set(gcf,'windowstyle','docked')
    legend('Original','Manipulated')
end






end

