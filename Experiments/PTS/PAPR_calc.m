function [ PAPR ] = PAPR_calc( x,L )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if isempty(L) % see Jiang& Wu figure 1
    L=4;
end

% Interpolation: see Jiang&Wu; real PAPR is obtained by zero padding in the
% frequency domain prior to IDFT.
% However, zero padding in the frequency domain is equivalent 
% to interpolating in the time domain. I prefer to interpolate in the time domain,
% so to enable comparison between original and manipulated signal
x_interp=interp(x,L); 


P_max=max(x_interp.*conj(x_interp));
P_avg=mean(x_interp.*conj(x_interp));

PAPR=db(P_max/P_avg,'power');


end

