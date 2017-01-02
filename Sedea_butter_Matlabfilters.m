%
% Second Order Butterworth filter coefficient Calculator Matlab Version
%
% Made by S Durbridge
%
% Last Edited: 02/01/2017
%
% Next Task: Add BPF and BSF
%

classdef Sedea_butter_Matlabfilters
    properties
        C
        A
        w0
        alpha
        Q
        fc
        fs
        gain
        a0
        a1
        a2
        b0
        b1
        b2
    end
    
    methods
        function obj = Sedea_butter_Matlabfilters(fc, fs, Q, gain)
            obj.fc = fc;
            obj.fs = fs;
            obj.Q = Q;
            obj.gain = gain;
        end
        function [coefs] = sedea_butter_lpf(obj)
            
            A = sqrt(10^(obj.gain/20));
            w0 = 2 * pi * obj.fc / obj.fs;
            C = 1/tan(w0);
            alpha = tan(w0) / (2 * obj.Q);
            
            a0 = 1/(1 + sqrt(2) * C + C^2);
            a1 = 2*a0;
            a2 = a0;
            b0 = A;
            b1 = 2 * a0 * (1 - C^2);
            b2 = a0 * (1 - sqrt(2) * C + C^2);
            
            
            b0 = b0 / a0;
            b1 =  b1 / a0;
            b2 =  b2 / a0;
            a1 =  a1 / a0;
            a2 =  a2 / a0;
            
            coefs = ([b0 b1 b2; a0 a1 a2]);
        end
        
        function  [coefs] = sedea_butter_hpf(obj)
            
            A = sqrt(10^(obj.gain/20));
            w0 = 2 * pi * obj.fc / obj.fs;
            C = tan(w0);
            alpha = tan(w0) / (2 * obj.Q);
            
            a0 = 1/(1 + sqrt(2) * C + C^2);
            a1 = -2*a0;
            a2 = a0;
            b0 = A;
            b1 = 2 * a0 * (C^2 - 1);
            b2 = a0 * (1 - sqrt(2) * C + C^2);
            
            
            b0 = b0 / a0;
            b1 =  b1 / a0;
            b2 =  b2 / a0;
            a1 =  a1 / a0;
            a2 =  a2 / a0;
            
            coefs = ([b0 b1 b2; a0 a1 a2]);
        end       
    end
end