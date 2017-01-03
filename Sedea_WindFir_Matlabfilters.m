%
% Fir filter coefficient Calculator Matlab Version
% Blackman Windowed
%
% Made by S Durbridge
%
% Last Edited: 03/01/2017
%
% Next Task: Add BPF and BSF and different windows
%

classdef Sedea_WindFir_Matlabfilters
    properties
      fs
      fc
      gain
      N
      df
      dfn
      fcn
      fc1
      fcn1
      hN
      n
      hd
      wd
      coefs
    end
    
    methods
        function obj = Sedea_WindFir_Matlabfilters(fc, fs, df, fc1)
            if((nargin > 0)&&(nargin < 4))
            obj.fc = fc;
            obj.fs = fs;
            obj.df = df;
            
            else
            obj.fc = fc;
            obj.fs = fs;
            obj.df = df;
            obj.fc1 = fc1;
            end   
        end
        function [coefs] = sedea_windfir_lpf(obj)
            
           dfn = obj.df/obj.fs;
           fcn = obj.fc/obj.fs;
           
           N = ceil(5.5/dfn);
           
           if (mod(N,2))
               hN = (N-1)/2;
               n = -hN:hN;
           else
               hN = N/2;
               n = -hN:hN-1;
           end
           
           hd = 2 * fcn * sin(2 * pi * fcn * n) ./ (n*2*pi*fcn);
           hd(hN +1) = 2 * fcn;
           
           wd = 0.42 + 0.5 * cos(2 * pi * n / (N-1)) + 0.08 * cos(4 * pi * n / (N-1));
           
           coefs = hd .* wd;
        end
        
        function [coefs] = sedea_windfir_hpf(obj)
            
           dfn = obj.df/obj.fs;
           fcn = obj.fc/obj.fs;
           
           N = ceil(5.5/dfn);
           
           if (mod(N,2))
               hN = (N-1)/2;
               n = -hN:hN;
           else
               hN = N/2;
               n = -hN:hN-1;
           end
           
           hd = -2 * fcn * sin(2 * pi * fcn * n) ./ (n*2*pi*fcn);
           hd(hN +1) = 1 - 2 * fcn;
           
           wd = 0.42 + 0.5 * cos(2 * pi * n / (N-1)) + 0.08 * cos(4 * pi * n / (N-1));
           
           coefs = hd .* wd;
        end 
        
        function [coefs] = sedea_windfir_bpf(obj)
            
           dfn = obj.df/obj.fs;
           fcn = obj.fc/obj.fs;
           fcn1 = obj.fc1/obj.fs;
           
           N = ceil(5.5/dfn);
           
           if (mod(N,2))
               hN = (N-1)/2;
               n = -hN:hN;
           else
               hN = N/2;
               n = -hN:hN-1;
           end
           
           hd = (2 * fcn1 * sin(2 * pi * fcn1 * n) ./ (n*2*pi*fcn1)) - (2 * fcn * sin(2 * pi * fcn * n) ./ (n*2*pi*fcn));
           hd(hN +1) = 2 * (fcn1 - fcn);
           
           wd = 0.42 + 0.5 * cos(2 * pi * n / (N-1)) + 0.08 * cos(4 * pi * n / (N-1));
           
           coefs = hd .* wd;
        end 
        function [coefs] = sedea_windfir_bsf(obj)
            
           dfn = obj.df/obj.fs;
           fcn = obj.fc/obj.fs;
           fcn1 = obj.fc1/obj.fs;
           
           N = ceil(5.5/dfn);
           
           if (mod(N,2))
               hN = (N-1)/2;
               n = -hN:hN;
           else
               hN = N/2;
               n = -hN:hN-1;
           end
           
           hd = (2 * fcn * sin(2 * pi * fcn * n) ./ (n*2*pi*fcn)) - (2 * fcn1 * sin(2 * pi * fcn1 * n) ./ (n*2*pi*fcn1));
           hd(hN +1) = 1 - 2 * (fcn - fcn1);
           
           wd = 0.42 + 0.5 * cos(2 * pi * n / (N-1)) + 0.08 * cos(4 * pi * n / (N-1));
           
           coefs = hd .* wd;
        end 
    end
end