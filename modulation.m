% function y = modulation(x,f_mod,f_sampling)
% Modulation of a signal
% Input : 
%   - signal x
%   - frequency of the modulation f_mod
%   - sampling frequency f_sampling
% Output :
%   - signal modulated y
function y = modulation(x,f_mod,f_sampling)

t = (0:length(x)-1) * (1/f_sampling);
g = exp(2*pi*j*f_mod*t);

y = x.*g;


