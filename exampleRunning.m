%Generate some filter coefs and do some plotting
clear all;
close all;
clc;
figure();
filc = Sedea_Rbj_Matlabfilters(10000, 48000, 2.0, 10.0);
somenums = sedea_rbjM_lpf(filc);
subplot(2,2,1);
zplane(somenums(1,:),somenums(2,:));
title('lpf coefs');
legend();
somenums1 = sedea_rbjM_hpf(filc);
subplot(2,2,2);
zplane(somenums1(1,:),somenums1(2,:));
title('hpf coefs');
legend();
somenums2 = sedea_rbjM_bpfcq(filc);
subplot(2,2,3);
zplane(somenums2(1,:),somenums1(2,:));
title('bpf coefs');
legend();
somenums3 = sedea_rbjM_bpfcg(filc);
somenums4 = sedea_rbjM_notch(filc);
subplot(2,2,4);
zplane(somenums4(1,:),somenums4(2,:));
title('notch coefs');
legend();
somenums5 = sedea_rbjM_apf(filc);
somenums6 = sedea_rbjM_pek(filc);
somenums7 = sedea_rbjM_ls(filc);
somenums8 = sedea_rbjM_hs(filc);

%generate a noise signal
Hpink = dsp.ColoredNoise(1,48000,1);
% Hsa = dsp.SpectrumAnalyzer('SampleRate',48000,'SpectrumType','Power density', ...
%     'OverlapPercent',50,'Window','Hamming', ...
%     'SpectralAverages',50,'PlotAsTwoSidedSpectrum',false, ...
%     'FrequencyScale','log','YLimits',[-50 20]);
pink = Hpink();
% Hsa(pink);
%Filter the noise
Hso = dsp.SpectrumAnalyzer('SampleRate',48000,'SpectrumType','Power density', ...
    'OverlapPercent',50,'Window','Hamming', ...
    'SpectralAverages',50,'PlotAsTwoSidedSpectrum',false, ...
    'FrequencyScale','log','YLimits',[-50 20]);

myout = Sedea_DFO_Matlabfilters(pink, somenums1(1,:), somenums1(2,:));
myout2 = Sedea_DFTWO_Matlabfilters(pink, somenums1(1,:), somenums1(2,:));
% pause(10.0);
myinfft = fft(pink);
myoutfft = fft(myout);
myoutfft2 = fft(myout2);

% Y = fft(X);
% Compute the two-sided spectrum P2. Then compute the single-sided spectrum P1 based on P2 and the even-valued signal length L.

myinfftP2 = abs(myinfft/length(myinfft));
myoutfftP2 = abs(myoutfft/length(myoutfft));
myoutfft2P2 = abs(myoutfft2/length(myoutfft2));
myinfftP1 = myinfftP2(1:size(myinfftP2)/2+1);
myoutfftP1 = myoutfftP2(1:size(myoutfftP2)/2+1);
myoutfft2P1 = myoutfft2P2(1:size(myoutfft2P2)/2+1);
myinfftP1(2:end-1) = 2*myinfftP1(2:end-1);
myoutfftP1(2:end-1) = 2*myoutfftP1(2:end-1);
myoutfft2P1(2:end-1) = 2*myoutfft2P1(2:end-1);
% Define the frequency domain f and plot the single-sided amplitude spectrum P1. The amplitudes are not exactly at 0.7 and 1, as expected, because of the added noise. On average, longer signals produce better frequency approximations.

f = 48000*(0:(length(myinfftP1)-1))/length(myinfftP1);
figure();
subplot(3,1,1);
plot(f,myinfftP1) ;
title('Single-Sided Amplitude Spectrum of Noise Signal');
xlabel('f (Hz)');
ylabel('|P1(f)|');
set(gca,'xscale','log');
subplot(3,1,2);
plot(f,myoutfftP1) ;
title('Single-Sided Amplitude Spectrum of Firect Form 1 Filtered Noise');
xlabel('f (Hz)');
ylabel('|P1(f)|');
set(gca,'xscale','log');
subplot(3,1,3);
plot(f,myoutfft2P1) ;
title('Single-Sided Amplitude Spectrum of Direct Form 2 Filtered Noise');
xlabel('f (Hz)');
ylabel('|P1(f)|');
set(gca,'xscale','log');
% Hso(myout);
% fvtool(somenums1(1,:), somenums1(2,:));

%plot the outputs
% figure();
figure();
afir = Sedea_WindFir_Matlabfilters(5000, 48000, 4000, 8000);
firlpf = sedea_windfir_bsf(afir);
freqz(firlpf, 66, 48000);
