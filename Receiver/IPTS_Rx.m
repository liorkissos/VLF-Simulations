function [ OFDM_matrix_PTS_Rx_f ] = IPTS_Rx( OFDM_matrix_Rx_f,OFDM_config)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% User defined params

Nguard_band_left=OFDM_config.Nguard_band_left;
Nguard_band_right=OFDM_config.Nguard_band_right;

P=size(OFDM_matrix_Rx_f,2); % number of OFDM symbols in the stream=columns in the OFDM matrix
N=size(OFDM_matrix_Rx_f,1)+Nguard_band_left+Nguard_band_right; % number of QAM symbols in the OFDM symbol=rows in the OFDM matrix


scrambling=OFDM_config.PTS.scrambling;
M=OFDM_config.PTS.M_PTS; % number of blocks the symbol is divided into

b_opt_mat=OFDM_config.PTS.b_opt_mat;


OFDM_matrix_Rx_w_GB_IPTS_f=zeros(N,P);

%%% need to convet to the exact form in which the PTS was applied in the Tx
%%% chain
OFDM_matrix_Rx_w_GB_PTS_f=[zeros(Nguard_band_left,P);OFDM_matrix_Rx_f;zeros(Nguard_band_right,P)]; %adding back the guard bands in order to correctly split the subcarriers
OFDM_matrix_Rx_w_GB_PTS_f=ifftshift(OFDM_matrix_Rx_w_GB_PTS_f,1);


for pp=1:P
    
    X_Rx=OFDM_matrix_Rx_w_GB_PTS_f(:,pp);
    
    b_opt=b_opt_mat(:,pp); % the PTS coefficients
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
    
end

OFDM_matrix_PTS_Rx_f=fftshift(OFDM_matrix_Rx_w_GB_IPTS_f,1); % bringing the symbol back to its initial arrangement at the function's input
OFDM_matrix_PTS_Rx_f=OFDM_matrix_PTS_Rx_f(Nguard_band_left+1:N-Nguard_band_right,:); % cutting out the guard bands



end

