%
% RBJ filter coefficient Calculator Matlab Version
%
% Made by S Durbridge
%
% Last Edited: 02/01/2017
%
% Next Task: Polish up as necesarry
%

classdef sedea_rbjM_filts
    
    methods

    function [coefs] = sedea_rbjM_lpf(fc, fs, Q, gain)
    
    A = sqrt(10^(gain/20));
    w0 = 2 * pi * fc / fs;
    alpha = sin(w0) / (2 * Q);
    
    a0 = A + alpha;
    b0 = (1 - cos(w0)) / 2;
    b1 = 1 - cos(w0);
    b2 = (1 - cos(w0)) / 2;
    a1 = -2 * cos(w0);
    a2 = 1 - alpha;
    
    b0 = b0 / a0;
    b1 =  b1 / a0;
    b2 =  b2 / a0;
    a1 =  a1 / a0;
    a2 =  a2 / a0;
    
    coefs = ([b0 b1 b2; a0 a1 a2]);
    end
    
    function  [coefs] = sedea_rbjM_hpf(fc, fs, Q, gain)
    
    
    A = sqrt(10^(gain/20));
    w0 = 2 * pi * fc / fs;
    alpha = sin(w0) / (2 * Q);
    
    a0 = A + alpha;
    b0 = (1 + cos(w0)) / 2;
    b1 = -(1 + cos(w0));
    b2 = (1 + cos(w0)) / 2;
    a1 = -2 * cos(w0);
    a2 = 1 - alpha;
    
    b0 = b0 / a0;
    b1 =  b1 / a0;
    b2 =  b2 / a0;
    a1 =  a1 / a0;
    a2 =  a2 / a0;
    
    coefs = ([b0 b1 b2; a0 a1 a2]);
    end
    
    function  [coefs] = sedea_rbjM_bpfcq(fc, fs, Q, gain)
    
    A = sqrt(10^(gain/20));
    w0 = 2 * pi * fc / fs;
    alpha = sin(w0) / (2 * Q);
    
    a0 = A + alpha;
    b0 = sin(w0)/2;
    b1 = 0.0;
    b2 = -sin(w0)/2;
    a1 = -2 * cos(w0);
    a2 = A - alpha;
    
    b0 = b0 / a0;
    b1 =  b1 / a0;
    b2 =  b2 / a0;
    a1 =  a1 / a0;
    a2 =  a2 / a0;
    
    coefs = ([b0 b1 b2; a0 a1 a2]);
    end
    
    function  [coefs] =  sedea_rbj_bpfcg(fc, fs, Q, gain)
    
    A = sqrt(10^(gain/20));
    w0 = 2 * pi * fc / fs;
    alpha = sin(w0) / (2 * Q);
    
    a0 = 1.0 + alpha;
    b0 = alpha;
    b1 = 0;
    b2 = -alpha;
    a1 = -2.0*cos(w0);
    a2 =  1.0 - alpha;
    
    b0 = b0 / a0;
    b1 =  b1 / a0;
    b2 =  b2 / a0;
    a1 =  a1 / a0;
    a2 =  a2 / a0;
    
    coefs = ([b0 b1 b2; a0 a1 a2]);
    end
    
    function  [coefs] = sedea_rbjM_notch(fc, fs, Q, gain)
    
    A = sqrt(10^(gain/20));
    w0 = 2 * pi * fc / fs;
    alpha = sin(w0) / (2 * Q);
    
    a0 = A + alpha;
    b0 = 1.0;
    b1 = -2.0*cos(w0);
    b2 = 1 + alpha;
    a1 = -2.0*cos(w0);
    a2 =  1.0 - alpha;
    
    b0 = b0 / a0;
    b1 =  b1 / a0;
    b2 =  b2 / a0;
    a1 =  a1 / a0;
    a2 =  a2 / a0;
    
    coefs = ([b0 b1 b2; a0 a1 a2]);
    end
    
    function  [coefs] =  sedea_rbjM_apf(fc, fs, Q, gain)
    
    A = sqrt(10^(gain/20));
    w0 = 2 * pi * fc / fs;
    alpha = sin(w0) / (2 * Q);
    
    b0 = 1.0 - alpha;
    b1 = -2.0*cos(w0);
    b2 = 1 + alpha;
    a0 = A + alpha;
    a1 = -2.0*cos(w0);
    a2 =  1.0 - alpha;
    
    b0 = b0 / a0;
    b1 =  b1 / a0;
    b2 =  b2 / a0;
    a1 =  a1 / a0;
    a2 =  a2 / a0;
    
    coefs = ([b0 b1 b2; a0 a1 a2]);
    end
    
    function  [coefs] = sedea_rbjM_pek(fc, fs, Q, gain)
    
    A = sqrt(10^(gain/20));
    w0 = 2 * pi * fc / fs;
    alpha = sin(w0) / (2 * Q);
    
    b0 = 1.0 + alpha * A;
    b1 = -2.0*cos(w0);
    b2 = 1 - alpha * A;
    a0 = 1 + alpha * A;
    a1 = -2.0*cos(w0);
    a2 =  1.0 - alpha / A;
    
    b0 = b0 / a0;
    b1 =  b1 / a0;
    b2 =  b2 / a0;
    a1 =  a1 / a0;
    a2 =  a2 / a0;
    
    coefs = ([b0 b1 b2; a0 a1 a2]);
    end
    
    function  [coefs] = sedea_rbjM_ls(fc, fs, Q, gain)
    
    A = sqrt(10^(gain/20));
    w0 = 2 * pi * fc / fs;
    alpha = sin(w0) / (2 * Q);
    
    b0 =      A * ((A+1) - (A-1)*cos(w0) + 2 * sqrt(A) * alpha);
    b1 =  2 * A * ((A-1) - (A+1)*cos(w0));
    b2 =      A * ((A+1) - (A-1)*cos(w0) - 2 * sqrt(A) * alpha);
    a0 =          ((A+1) + (A-1)*cos(w0) + 2 * sqrt(A) * alpha);
    a1 = -2 *     ((A-1) + (A+1)*cos(w0));
    a2 =          ((A+1) + (A-1)*cos(w0) - 2 * sqrt(A) * alpha);
    
    b0 = b0 / a0;
    b1 =  b1 / a0;
    b2 =  b2 / a0;
    a1 =  a1 / a0;
    a2 =  a2 / a0;
    
    coefs = ([b0 b1 b2; a0 a1 a2]);
    end
    
    function  [coefs] = sedea_rbj_hs(fc, fs, Q, gain)
    
    A = sqrt(10^(gain/20));
    w0 = 2 * pi * fc / fs;
    alpha = sin(w0) / (2 * Q);
    
    b0 =      A *  ((A+1) - (A-1)*cos(w0) + 2 * sqrt(A) * alpha);
    b1 = -2 * A *  ((A-1) - (A+1)*cos(w0));
    b2 =      A *  ((A+1) - (A-1)*cos(w0) - 2 * sqrt(A) * alpha);
    a0 =           ((A+1) + (A-1)*cos(w0) + 2 * sqrt(A) * alpha);
    a1 =  2     *  ((A-1) + (A+1)*cos(w0));
    a2 =           ((A+1) + (A-1)*cos(w0) - 2 * sqrt(A) * alpha);
    
    b0 = b0 / a0;
    b1 =  b1 / a0;
    b2 =  b2 / a0;
    a1 =  a1 / a0;
    a2 =  a2 / a0;
    
    coefs = ([b0 b1 b2; a0 a1 a2]);
    end
end
% somenums[] = sedea_rbjM_lpf(10000, 48000, 2.0, -1.0);
% somenums1 = sedea_rbjM_hpf(10000, 48000, 3.0, -1.0);
% somenums2 = sedea_rbjM_bpfcq(10000, 48000, 3.0, -1.0);
% somenums3 = sedea_rbjM_bpfcg(5000, 48000, 3.0, -1.0);
% somenums4 = sedea_rbjM_notch(5000, 48000, 3.0, -1.0);
% somenums5 = sedea_rbjM_apf(5000, 48000, 10, 1.0);
% somenums6 = sedea_rbjM_pek(5000, 48000, 10, 1.0);
% somenums7 = sedea_rbjM_ls(5000, 48000, 1.0, -1.0);
% somenums8 = sedea_rbjM_hs(5000, 48000, 1.0, -1.0);
% % #b = signal.firwin(80, 0.5, window=('kaiser', 8))
%
% filterAnal(somenums[0], somenums[1], 48000);
% filterAnal(somenums1[0], somenums1[1], 48000);
% filterAnal(somenums2[0], somenums2[1], 48000)#;
% filterAnal(somenums3[0], somenums3[1], 48000);
% % filterAnal(somenums4[0], somenums4[1], 48000);
% % filterAnal(somenums5[0], somenums5[1], 48000);
% % filterAnal(somenums6[0], somenums6[1], 48000);
% % filterAnal(somenums7[0], somenums7[1], 48000);
% filterAnal(somenums8[0], somenums8[1], 48000)
% '''
% w, h = signal.freqz(b, a)
%
%
% fig = plt.figure()
% plt.title('Digital filter frequency response')
% ax1 = fig.add_subplot(111)
%
% plt.plot(w, 20 * np.log10(abs(h)), 'b')
% plt.ylabel('Amplitude [dB]', color='b')
% plt.xlabel('Frequency [rad/sample]')
%
% ax2 = ax1.twinx()
% angles = np.unwrap(np.angle(h))
% plt.plot(w, angles, 'g')
% plt.ylabel('Angle (radians)', color='g')
% plt.grid()
% plt.axis('tight')
% plt.show()
% '''