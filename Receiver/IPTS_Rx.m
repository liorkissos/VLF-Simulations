function [ OFDM_matrix_PTS_Rx_f ] = IPTS_Rx( OFDM_matrix_Rx_f,OFDM_config)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% User defined params

%PTS_algorithm=OFDM_config.PTS.PTS_algorithm;
scrambling=OFDM_config.PTS.scrambling;

M=OFDM_config.PTS.M_PTS; % number of blocks the symbol is divided into
%W=OFDM_config.PTS.W_PTS; % number of possible phases
%L=OFDM_config.PTS.L_PTS; % interpolation factor. see Jiang& Wu equation 5

%P_data=OFDM_config.PTS.P_data; %required power from
%N_preamble_CE=OFDM_config.N_preamble_CE;

%scrambling=OFDM_config.PTS.scrambling;
b_opt_mat=OFDM_config.PTS.b_opt_mat;

P=size(OFDM_matrix_Rx_f,2); % number of OFDM symbols in the stream=columns in the OFDM matrix

Nguard_band_left=OFDM_config.Nguard_band_left;
Nguard_band_right=OFDM_config.Nguard_band_right;

N=size(OFDM_matrix_Rx_f,1)+Nguard_band_left+Nguard_band_right; % number of QAM symbols in the OFDM symbol=rows in the OFDM matrix

OFDM_matrix_PTS_Rx_f=zeros(size(OFDM_matrix_Rx_f,1),P);

%%% need to convet to the exact form in which the PTS was applied 
OFDM_matrix_Rx_w_GB_f=[zeros(Nguard_band_left,P);OFDM_matrix_Rx_f;zeros(Nguard_band_right,P)]; %adding back the guard bands in order to correctly split the subcarriers 
OFDM_matrix_Rx_w_GB_f=ifftshift(OFDM_matrix_Rx_w_GB_f,1);
%%%

for pp=1:P
    
    
    %X_Rx=[zeros(Nguard_band_left,1);OFDM_matrix_Rx_f(:,pp);zeros(Nguard_band_right,1)]; %adding back the guard bands in order to correctly split the subcarriers 
    X_Rx=OFDM_matrix_Rx_w_GB_f(:,pp);

    
    b_opt=b_opt_mat(:,pp);
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
    
    
    X_IPTS_Rx=sum(X_Rx_mat,2);
    OFDM_matrix_Rx_w_GB_IPTS_f(:,pp)=X_IPTS_Rx;
   % OFDM_matrix_PTS_Rx_f(:,pp)=X_IPTS_Rx(Nguard_band_left+1:N-Nguard_band_right);
    
end

OFDM_matrix_PTS_Rx_f=fftshift(OFDM_matrix_Rx_w_GB_IPTS_f,1);
OFDM_matrix_PTS_Rx_f=OFDM_matrix_PTS_Rx_f(Nguard_band_left+1:N-Nguard_band_right,:);



end

