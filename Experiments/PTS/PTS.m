
clc
clear
close all


% %%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% Complex
% %%%%%%%%%%%%%%%%%%%%%%%%%
% 
% M=4;
% N=2^6;
% 
% W=4;
% 
% Modulation=4;
% 
% 
% %% Data Generation
% 
% X_orig=randi(Modulation,N,1);
% 
% X_orig_shift=ifftshift(X_orig);
% 
% 
% %% Phase optimization
% 
% phase=randi(W,M)-1;
% b=exp(j*2*pi*phase/W);
% 
% % b=ones(M,1);
% % 
% % for kk=2:M
% %     
% %     for ll=0:W-1
% %         
% %     end
% %     
% % end
% 
% %% Tx
% 
% 
% for kk=1:M
%     ind_vec=kk:M:N;
%     X_subs=X_orig_shift(ind_vec);
%     X_orig_mat(:,kk)=b(kk)*upsample(X_subs,M,kk-1);
%     
% end
% 
% X_Tx_mat=X_orig_mat;
% x_Tx_mat=ifft(X_Tx_mat);
% 
% x_Tx=sum(x_Tx_mat,2);
% 
% 
% %% PAPR measurement
% 
% x_Tx_orig=ifft(X_orig_shift);
% P_max_orig=max(abs(x_Tx_orig));
% P_avg_orig=mean(abs(x_Tx_orig));
% PAPR_orig=10*log10(P_max_orig/P_avg_orig)
% 
% P_max_manip=max(abs(x_Tx));
% P_avg_manip=mean(abs(x_Tx));
% PAPR_manip=10*log10(P_max_manip/P_avg_manip)
% 
% %% Rx
% 
% X_Rx=fft(x_Tx);
% 
% for kk=1:M
%     ind_vec=kk:M:N;
%     X_Rx_subs=X_Rx(ind_vec);
%     X_Rx_mat(:,kk)=(1./b(kk))*upsample(X_Rx_subs,M,kk-1);
%     
% end
% 
% X_Rx=sum(X_Rx_mat,2);
% 
% X_Rx=fftshift(X_Rx);
% 
% error=max(abs(X_orig-X_Rx))
% figure
% set(gcf,'windowstyle','docked')
% shg
% stem(real(X_orig));hold on;stem(real(X_Rx))
% title(['Error= ',num2str(error),''])
% legend('original','reconstructed')
% grid minor



%%%%%%%%%%%%%%%%%%%%%%%%
%%% Final
%%%%%%%%%%%%%%%%%%%5


M=4;
N=2^6;

W=4;

Modulation=4;


%% Data Generation

X_orig=randi(Modulation,N,1);

X_orig_shift=ifftshift(X_orig);


%% Phase optimization

phase=randi(W,M)-1;
b=exp(j*2*pi*phase/W);

% b=ones(M,1);
% 
% for kk=2:M
%     
%     for ll=0:W-1
%         
%     end
%     
% end

%% Tx


for kk=1:M
    ind_vec=kk:M:N;
    X_subs=X_orig_shift(ind_vec);
    X_orig_mat(:,kk)=b(kk)*upsample(X_subs,M,kk-1);
    
end

X_Tx_mat=X_orig_mat;
x_Tx_mat=ifft(X_Tx_mat);

x_Tx=sum(x_Tx_mat,2);


%% PAPR measurement

x_Tx_orig=ifft(X_orig_shift);
P_max_orig=max(abs(x_Tx_orig));
P_avg_orig=mean(abs(x_Tx_orig));
PAPR_orig=10*log10(P_max_orig/P_avg_orig)

P_max_manip=max(abs(x_Tx));
P_avg_manip=mean(abs(x_Tx));
PAPR_manip=10*log10(P_max_manip/P_avg_manip)

%% Rx

X_Rx=fft(x_Tx);

for kk=1:M
    ind_vec=kk:M:N;
    X_Rx_subs=X_Rx(ind_vec);
    X_Rx_mat(:,kk)=(1./b(kk))*upsample(X_Rx_subs,M,kk-1);
    
end

X_Rx=sum(X_Rx_mat,2);

X_Rx=fftshift(X_Rx);

error=max(abs(X_orig-X_Rx))
figure
set(gcf,'windowstyle','docked')
shg
stem(real(X_orig));hold on;stem(real(X_Rx))
title(['Error= ',num2str(error),''])
legend('original','reconstructed')
grid minor


