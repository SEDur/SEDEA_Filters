%
% RBJ filter coefficient Calculator Matlab Version
%
% Made by S Durbridge
%
% Last Edited: 02/01/2017
%
% Next Task: Polish up as necesarry
%

classdef Sedea_Rbj_Matlabfilters
    properties
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
        function obj = Sedea_Rbj_Matlabfilters(fc, fs, Q, gain)
            obj.fc = fc;
            obj.fs = fs;
            obj.Q = Q;
            obj.gain = gain;
        end
        function [coefs] = sedea_rbjM_lpf(obj)
            
            A = sqrt(10^(obj.gain/20));
            w0 = 2 * pi * obj.fc / obj.fs;
            alpha = sin(w0) / (2 * obj.Q);
            
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
        
        function  [coefs] = sedea_rbjM_hpf(obj)
            
            
            A = sqrt(10^(obj.gain/20));
            w0 = 2 * pi * obj.fc / obj.fs;
            alpha = sin(w0) / (2 * obj.Q);
            
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
        
        function  [coefs] = sedea_rbjM_bpfcq(obj)
            
            A = sqrt(10^(obj.gain/20));
            w0 = 2 * pi * obj.fc / obj.fs;
            alpha = sin(w0) / (2 * obj.Q);
            
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
        
        function  [coefs] =  sedea_rbjM_bpfcg(obj)
            
            A = sqrt(10^(obj.gain/20));
            w0 = 2 * pi * obj.fc / obj.fs;
            alpha = sin(w0) / (2 * obj.Q);
            
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
        
        function  [coefs] = sedea_rbjM_notch(obj)
            
            A = sqrt(10^(obj.gain/20));
            w0 = 2 * pi * obj.fc / obj.fs;
            alpha = sin(w0) / (2 * obj.Q);
            
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
        
        function  [coefs] =  sedea_rbjM_apf(obj)
            
            A = sqrt(10^(obj.gain/20));
            w0 = 2 * pi * obj.fc / obj.fs;
            alpha = sin(w0) / (2 * obj.Q);
            
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
        
        function  [coefs] = sedea_rbjM_pek(obj)
            
            A = sqrt(10^(obj.gain/20));
            w0 = 2 * pi * obj.fc / obj.fs;
            alpha = sin(w0) / (2 * obj.Q);
            
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
        
        function  [coefs] = sedea_rbjM_ls(obj)
            
            A = sqrt(10^(obj.gain/20));
            w0 = 2 * pi * obj.fc / obj.fs;
            alpha = sin(w0) / (2 * obj.Q);
            
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
        
        function  [coefs] = sedea_rbjM_hs(obj)
            
            A = sqrt(10^(obj.gain/20));
            w0 = 2 * pi * obj.fc / obj.fs;
            alpha = sin(w0) / (2 * obj.Q);
            
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
end