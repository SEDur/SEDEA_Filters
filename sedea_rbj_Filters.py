# -*- coding: utf-8 -*-
"""
RBJ filter coefficient Calculator

Made by S Durbridge

Last Edited: 08/12/2016

Next Task: Write the framework for the coefficient calculator

"""
import numpy as np
from scipy import signal
import matplotlib.pyplot as plt


def sedea_rbj_lpf(fc, fs, Q, gain):
    
    b0 = 1.0
    A = np.sqrt(10**(gain/20))
    w0 = 2 * np.pi * fc / fs
    alpha = np.sin(w0) / (2 * Q)
    
    a0 = A + alpha
    b0 = (1 - np.cos(w0)) / 2
    b1 = 1 - np.cos(w0)
    b2 = (1 - np.cos(w0)) / 2
    a1 = -2 * np.cos(w0)
    a2 = 1 - alpha
    
    b0 /=  a0
    b1 /=  a0
    b2 /=  a0
    a1 /=  a0
    a2 /=  a0
        
    
    coefs = np.array([b0, b1, b2, a0, a1, a2]);    
    return (coefs)
    
def sedea_rbj_hpf(fc, fs, Q, gain):
    
    b0 = 1.0
    A = np.sqrt(10**(gain/20))
    w0 = 2 * np.pi * fc / fs
    alpha = np.sin(w0) / (2 * Q)
    
    a0 = A + alpha
    b0 = (1 + np.cos(w0)) / 2
    b1 = -(1 + np.cos(w0))
    b2 = (1 + np.cos(w0)) / 2
    a1 = -2 * np.cos(w0)
    a2 = 1 - alpha
       
    b0 /=  a0
    b1 /=  a0
    b2 /=  a0
    a1 /=  a0
    a2 /=  a0
    
    coefs = np.array([b0, b1, b2, a0, a1, a2]);    
    return(coefs)

somenums = sedea_rbj_lpf(500, 48000, 3, 90)
    
#b = signal.firwin(80, 0.5, window=('kaiser', 8))
b = np.array([somenums[0], somenums[1], somenums[2]])
a = np.array([somenums[3], somenums[4], somenums[5]])
w, h = signal.freqz(b, a)

fig = plt.figure()
plt.title('Digital filter frequency response')
ax1 = fig.add_subplot(111)

plt.plot(w, 20 * np.log10(abs(h)), 'b')
plt.ylabel('Amplitude [dB]', color='b')
plt.xlabel('Frequency [rad/sample]')

ax2 = ax1.twinx()
angles = np.unwrap(np.angle(h))
plt.plot(w, angles, 'g')
plt.ylabel('Angle (radians)', color='g')
plt.grid()
plt.axis('tight')
plt.show()
