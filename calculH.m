function H = calculH(x,ordre,memoire,f0,fs)

%% Build the matrix
T = toeplitz(x');
H = repmat(T(:,1:memoire),1,ordre)';

%% Modification for harmonic and subharmonic components
for kO = 2:ordre
    for kSignal=(kO-1)*memoire+1:kO*memoire
        H(kSignal,:) = real(modulation(hilbert(H(kSignal,:)),(kO-1)*f0,fs));
    end
end

