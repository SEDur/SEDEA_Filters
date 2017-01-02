%
% Second Order Butterworth filter coefficient Calculator Matlab Version
%
% Made by S Durbridge
%
% Last Edited: 02/01/2017
%
% Next Task: Add BPF and BSF
%

classdef Sedea_LRF_Matlabfilters
    properties
        Theta
        Kappa
        Omega
        Delta
        A
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
        function obj = Sedea_LRF_Matlabfilters(fc, fs, Q, gain)
            obj.fc = fc;
            obj.fs = fs;
            obj.Q = Q;
            obj.gain = gain;
        end
        function [coefs] = sedea_LRF_lpf(obj)
            
            A = sqrt(10^(obj.gain/20));
            Theta = pi * obj.fc / obj.fs;
            Omega = pi * obj.fc;
            Kappa = Omega / tan(Theta);
            Delta = Kappa^2 + Omega^2 + 2*Kappa*Omega;
            
            a0 = Omega^2 / Delta;
            a1 = 2 * (Omega^2 / Delta);
            a2 = Omega^2/Delta;
            b0 = A;
            b1 = (-2*Kappa^2 + 2*Omega^2)/Delta;
            b2 = (-2*Kappa*Omega + Kappa^2 + Omega^2)/Delta;
            
            
            b0 = b0 / a0;
            b1 =  b1 / a0;
            b2 =  b2 / a0;
            a1 =  a1 / a0;
            a2 =  a2 / a0;
            
            coefs = ([b0 b1 b2; a0 a1 a2]);
        end
        
        function  [coefs] = sedea_LRF_hpf(obj)
            
            A = sqrt(10^(obj.gain/20));
            Theta = pi * obj.fc / obj.fs;
            Omega = pi * obj.fc;
            Kappa = Omega / tan(Theta);
            Delta = Kappa^2 + Omega^2 + 2*Kappa*Omega;
            
            a0 = Kappa^2 / Delta;
            a1 = 2 * (Kappa^2 / Delta);
            a2 = Kappa^2/Delta;
            b0 = A;
            b1 = (-2*Kappa^2 + 2*Omega^2)/Delta;
            b2 = (-2*Kappa*Omega + Kappa^2 + Omega^2)/Delta;
            
            
            b0 = b0 / a0;
            b1 =  b1 / a0;
            b2 =  b2 / a0;
            a1 =  a1 / a0;
            a2 =  a2 / a0;
            
            coefs = ([b0 b1 b2; a0 a1 a2]);
        end       
    end
end