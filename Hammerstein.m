% function [z_model,z_decomp,T] = Hammerstein(x,z,ordre,memoire,f0,fs)
% Hammerstein reconstruction
% Input : 
%   - input signal x
%   - output signal z
%   - order of the Hammerstein structure ordre
%   - memory of the Hammertein structure memoire
%   - frequency of the subharmonic f0
%   - sampling frequency fs
% Output :
%   - Model of the output signal z_model
%   - Signals for each output components z_decomp

function [z_model,z_decomp,T] = Hammerstein(x,z,ordre,memoire,f0,fs)

%% Build the input matrix
H = calculH(x,ordre,memoire,f0,fs);

%%
Y = z';
R = H*H';
P = H*Y;

%%
[U,S,V] = svd(R);
D       = diag(S);
Ds      = sort(D,'descend');
index   = find(Ds>=0.01);
k       = index(end);
R_hat   = V(:,1:k)*pinv(S(1:k,1:k))*U(:,1:k)'; 
T       = R_hat*P;

%% Model of the output signal
z_model = (H' * T);

%% Signals for each output components
z_decomp = zeros(length(z), ordre);
for k = 1:ordre
    z_decomp(:, k) = H((k-1)*memoire+1:k*memoire,:)' * T((k-1)*memoire+1:k*memoire);
end
