function [ Dmin] = Aliasing_check( Frec,Fs,F_if,security )
%ALIASING_CHECK Summary of this function goes here
%   Detailed explanation goes here



Frec=Frec/1e3;
F_if=F_if/1e3;
Fs=Fs/1e3;
security=security/1e3;


a=[-3*Frec:Frec:-1*Frec,1*Frec:Frec:3*Frec] ; % only first 3 replicas are significant (saw on spectrum)
A=repmat(a,[31 1]);

B1=A+F_if; % D/A replicas: right hand side
B2=A-F_if; % D/A replicas: left hand side 

%C=repmat([(-15*Fs:Fs:-1*Fs),(1*Fs:Fs:15*Fs)]',[1 10]); % A/D's aliasing translation
C=repmat([(-15*Fs:Fs:15*Fs)]',[1 length(a)]); % A/D's aliasing translation: up to 15 Nyquist zones


% D/A's replicas aliases 
D1=C+B1;
D2=C+B2;

% D/A's replicas aliases into 1st nyquist zone: center frequencies
%[replica_1_row,replica_1_column]=find(D1>-(Fs/2) & D1<(Fs/2));

E1=D1(D1>-(Fs/2) & D1<(Fs/2));
E2=D2(D2>-(Fs/2) & D2<(Fs/2));

%[replica_1_row,replica_1_column]=find(max(abs(F_if-abs(E1))<security)


if max(abs(F_if-abs(E1))<security  | abs(F_if-abs(E2))<security )
    warning('D/A s replicas are aliased onto signal')
end
    
Dmin=min(min(abs(F_if-abs(E1))),min(abs(F_if-abs(E2)))); % minimum distance between F_if and aliaed replica's center



end

