
clc
clear
close all

%%% Simple
%%%%%%%%%%%%%%%%%%%5

% N=2^4;
% 
% M=N;
% 
% X_orig=1:1:N;
% 
% X=ifftshift(X_orig);
% x=ifft(X,M);
% 
% X_rec=fft(x);
% X_rec=fftshift(X_rec);
% 
% error=max(abs(X_orig-X_rec))
% stem(X_orig);hold on;stem(X_rec)
% legend('original','reconstructed')


%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Moderate
%%%%%%%%%%%%%%%%%%%%%%%%%

% K=2;
% N=2^4;
% 
% M=K*N;
% 
% 
% %% Tx
% 
% X_orig=1:1:N;
% X_orig=X_orig.';
% 
% X_orig_shift=ifftshift(X_orig);
% 
% for kk=1:K
%     ind_vec=kk:K:N;
%     X_subs=X_orig_shift(ind_vec);
%     X_orig_mat(:,kk)=upsample(X_subs,K,kk-1);
%     
% end
% 
% X_Tx_mat=X_orig_mat;
% x_Tx_mat=ifft(X_Tx_mat);
% 
% x_Tx=sum(x_Tx_mat,2);
% 
% %% Rx
% 
% X_Rx=fft(x_Tx);
% X_Rx=fftshift(X_Rx);
% 
% error=max(abs(X_orig-X_Rx))
% stem(X_orig);hold on;stem(X_Rx)
% legend('original','reconstructed')
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Complex
%%%%%%%%%%%%%%%%%%%%%%%%%

M=2;
N=2^3;

Modulation=4;



%% Tx

b=rand(M,1);

%X_orig=1:1:N;
%X_orig=X_orig.';

X_orig=randi(Modulation,N,1);

X_orig_shift=ifftshift(X_orig);

for kk=1:M
    ind_vec=kk:M:N;
    X_subs=X_orig_shift(ind_vec);
    X_orig_mat(:,kk)=b(kk)*upsample(X_subs,M,kk-1);
    
end

X_Tx_mat=X_orig_mat;
x_Tx_mat=ifft(X_Tx_mat);

x_Tx=sum(x_Tx_mat,2);

%% Rx

X_Rx=fft(x_Tx);

for kk=1:M
    ind_vec=kk:M:N;
    X_Rx_subs=X_Rx(ind_vec);
    X_Rx_mat(:,kk)=(1./b(kk))*upsample(X_Rx_subs,M,kk-1);
    
end

%x_Rx_mat=fft(X_Rx_mat);
X_Rx=sum(X_Rx_mat,2);

X_Rx=fftshift(X_Rx);

error=max(abs(X_orig-X_Rx))
stem(X_orig);hold on;stem(X_Rx)
legend('original','reconstructed')
