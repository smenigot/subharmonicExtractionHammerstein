% function y = modulation(x,f_mod,f_sampling)
% demodulation of a signal
% Input : 
%   - signal x
%   - frequency of the demodulation f_mod
%   - sampling frequency f_sampling
% Output :
%   - signal demodulated y

function y = demodulation(x,f_mod,f_sampling)

t = (0:length(x)-1) * (1/f_sampling);
g = conj(exp(2*pi*j*f_mod*t));

y = x.*g;


