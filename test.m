close all
clear
clc
%%
N  = 512;  % number of samples
fe = 100e6; % frequency sampling

t = (0:N-1) * 1/fe; % time

%% input signal 
f0 = 5e6; % transmit frequency
x  = cos(2*pi*f0.*t); 

%% output signal
y  = rand*cos(2*pi*f0.*t) + ... % fundamental
    rand*cos(2*pi*2*f0.*t) + ... % second harmonic
    rand*cos(2*pi*3*f0.*t) + ... % third harmonic 
    rand*cos(2*pi*f0/2.*t) + ... % subharmonic
    rand*cos(2*pi*3*f0/2.*t); % ultraharmonic

%% Modified Input and Output
x_demod = real(demodulation(hilbert(x),f0/2,fe));
y_demod = real(modulation(hilbert(y),f0/2,fe));

%% Hammerstein with Modified Input
ordre   = 2;            % order
memoire = length(y);    % memory

[y_modelX,y_decompX] = Hammerstein(x_demod, y, ordre, memoire, f0/2, fe);


%% Hammerstein with Modified Output
ordre   = 2;            % order
memoire = length(y);    % memory

[y_modelY,y_decompY] = Hammerstein(x, y_demod, ordre, memoire,f0/2,fe);

% Demodulation output to recover the original signal
y_modelY = real(demodulation(hilbert(y_modelY'),f0/2,fe));
for k=1:ordre
    y_decompY(:,k) = real(demodulation(hilbert(y_decompY(:,k)'),f0/2,fe))';
end

%% Standard Filtering for Comparison
wb1=2*4e6/fe;
wb2=2*3e6/fe;
wh1=2*6e6/fe;
wh2=2*7e6/fe;
Wp = [wb1 wh1];
Ws = [wb2 wh2];
[n2,Wn2] = buttord(Wp,Ws,3,10);

[B2,A2]=butter(n2,Wn2);

y_filt_sub=filtfilt(B2,A2,y);

%% Figure of Spectra
figure,
N = 2*8192;
freq = (0:N-1)/N * fe;
hl = subplot(3,1,1);
plot(freq/1e6,20*log10(abs(fft(x,N))/max(abs(fft(x,N)))),...
    freq/1e6,20*log10(abs(fft(y,N))/max(abs(fft(y,N)))),...
    'LineWidth',2)
set(hl,'FontWeight','b','FontSize',14)
xlabel('Frequency (MHz)')
ylabel('Amplitude (dB)')
xlim([0 16])
ylim([-25 0])
legend('Excitation','Microbubble','Location','NorthWest')
box on, grid on
title('(a)')

hl = subplot(3,1,2);
plot(freq/1e6,20*log10(abs(fft(x_demod,N))/max(abs(fft(x_demod,N)))),'--',...
    freq/1e6,20*log10(abs(fft(y,N))/max(abs(fft(y,N)))),...    
    freq/1e6,20*log10(abs(fft(y_modelY,N))/max(abs(fft(y_modelX,N)))),'-.',...
    'LineWidth',2)
set(hl,'FontWeight','b','FontSize',14)
xlabel('Frequency (MHz)')
ylabel('Amplitude (dB)')
xlim([0 16])
ylim([-25 0])
hl = legend('Input 1','Output 1','Model','Location','NorthWest')
box on, grid on
title('(b)')

hl = subplot(3,1,3);
plot(freq/1e6,20*log10(abs(fft(x,N))/max(abs(fft(x,N)))),...
    freq/1e6,20*log10(abs(fft(y_demod,N))/max(abs(fft(y_demod,N)))),'--',...
    freq/1e6,20*log10(abs(fft(y_modelY,N))/max(abs(fft(y_modelY,N)))),'-.',...
    'LineWidth',2)
set(hl,'FontWeight','b','FontSize',14)
xlabel('Frequency (MHz)')
ylabel('Amplitude (dB)')
xlim([0 16])
ylim([-25 0])
legend('Input 2','Output 2','Model after demodulation','Location','NorthWest')
box on, grid on
title('(c)')


%% Figure of Time Representation
figure,
hl = subplot(3,1,1);
plot(t/1e-6,y,...
    t/1e-6,y_modelX','--',...
    t/1e-6,y_modelY,'-.',...
    t/1e-6,y-y_modelX',...    
    t/1e-6,y-y_modelY,...
    'LineWidth',2)
set(hl,'FontWeight','b','FontSize',14)
title('(a)')
legend('Microbubble','Model 1','Model 2','error 1','error 2')
box on, grid on
xlabel('Time (\mus)')
ylabel('Amplitude (u. a.)')

hl = subplot(312);
plot(t/1e-6,y_decompX(:,1),'--',...
    t/1e-6,y_decompY(:,1),'-.',...
    t/1e-6,y_filt_sub,...
    'LineWidth',2)
set(hl,'FontWeight','b','FontSize',14)
title('(b)')
legend('Subharmonics from model 1','Subharmonics from model 2','Subharmonics from standard filtering','Location','SouthWest')
box on, grid on
xlabel('Time (\mus)')
ylabel('Amplitude (u. a.)')

N = 2*8192;
freq = (0:N-1)/N * fe;
hl = subplot(313);
plot(freq/1e6,20*log10(abs(fft(y_filt_sub,N))/max(abs(fft(y,N)))),...
    freq/1e6,20*log10(abs(fft(y_decompX(:,1),N))/max(abs(fft(y,N)))),'--',...
    freq/1e6,20*log10(abs(fft(y_decompY(:,1),N))/max(abs(fft(y,N)))),'-.',...
    freq/1e6,20*log10(abs(fft(y,N))/max(abs(fft(y,N)))),...
    'LineWidth',2)
set(hl,'FontWeight','b','FontSize',14)
xlim([0 12])
ylim([-25 0])
title('(c)')
legend('Subharmonics from model 1','Subharmonics from model 2','Subharmonics from standard filtering','Microbubble','Location','NorthWest')
box on, grid on
xlabel('Time (\mus)')
ylabel('Amplitude (dB)')

