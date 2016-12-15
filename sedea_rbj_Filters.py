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
        
    
    coefs = np.array([np.array([b0, b1, b2]), np.array([a0, a1, a2])]);    
    return (coefs)
    
def sedea_rbj_hpf(fc, fs, Q, gain):
    
    
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
    
    coefs = np.array([np.array([b0, b1, b2]), np.array([a0, a1, a2])]);   
    return(coefs)
    
def sedea_rbj_bpfcq(fc, fs, Q, gain):
    
    A = np.sqrt(10**(gain/20))
    w0 = 2 * np.pi * fc / fs
    alpha = np.sin(w0) / (2 * Q)
    
    a0 = A + alpha
    b0 = np.sin(w0)/2
    b1 = 0.0
    b2 = -np.sin(w0)/2
    a1 = -2 * np.cos(w0)
    a2 = A - alpha

    b0 /=  a0
    b1 /=  a0
    b2 /=  a0
    a1 /=  a0
    a2 /=  a0
    
    coefs = np.array([np.array([b0, b1, b2]), np.array([a0, a1, a2])]);    
    return(coefs)
    
def sedea_rbj_bpfcg(fc, fs, Q, gain):
    
    A = np.sqrt(10**(gain/20))
    w0 = 2 * np.pi * fc / fs
    alpha = np.sin(w0) / (2 * Q)
    
    a0 = 1.0 + alpha
    b0 = alpha
    b1 = 0
    b2 = -alpha
    a1 = -2.0*np.cos(w0)
    a2 =  1.0 - alpha

    b0 /=  a0
    b1 /=  a0
    b2 /=  a0
    a1 /=  a0
    a2 /=  a0
    
    coefs = np.array([np.array([b0, b1, b2]), np.array([a0, a1, a2])]);    
    return(coefs)
    
def sedea_rbj_notch(fc, fs, Q, gain):
    
    A = np.sqrt(10**(gain/20))
    w0 = 2 * np.pi * fc / fs
    alpha = np.sin(w0) / (2 * Q)
    
    a0 = A + alpha
    b0 = 1.0
    b1 = -2.0*np.cos(w0)
    b2 = 1 + alpha
    a1 = -2.0*np.cos(w0)
    a2 =  1.0 - alpha

    b0 /=  a0
    b1 /=  a0
    b2 /=  a0
    a1 /=  a0
    a2 /=  a0
    
    coefs = np.array([np.array([b0, b1, b2]), np.array([a0, a1, a2])]);    
    return(coefs)
    
def sedea_rbj_apf(fc, fs, Q, gain):
    
    A = np.sqrt(10**(gain/20))
    w0 = 2 * np.pi * fc / fs
    alpha = np.sin(w0) / (2 * Q)
    
    b0 = 1.0 - alpha
    b1 = -2.0*np.cos(w0)
    b2 = 1 + alpha
    a0 = A + alpha
    a1 = -2.0*np.cos(w0)
    a2 =  1.0 - alpha

    b0 /=  a0
    b1 /=  a0
    b2 /=  a0
    a1 /=  a0
    a2 /=  a0
    
    coefs = np.array([np.array([b0, b1, b2]), np.array([a0, a1, a2])]);    
    return(coefs)
    
def sedea_rbj_pek(fc, fs, Q, gain):
    
    A = np.sqrt(10**(gain/20))
    w0 = 2 * np.pi * fc / fs
    alpha = np.sin(w0) / (2 * Q)
    
    b0 = 1.0 + alpha * A
    b1 = -2.0*np.cos(w0)
    b2 = 1 - alpha * A
    a0 = 1 + alpha * A
    a1 = -2.0*np.cos(w0)
    a2 =  1.0 - alpha / A

    b0 /=  a0
    b1 /=  a0
    b2 /=  a0
    a1 /=  a0
    a2 /=  a0
    
    coefs = np.array([np.array([b0, b1, b2]), np.array([a0, a1, a2])]);    
    return(coefs)
    
def sedea_rbj_ls(fc, fs, Q, gain):
    
    A = np.sqrt(10**(gain/20))
    w0 = 2 * np.pi * fc / fs
    alpha = np.sin(w0) / (2 * Q)
    
    b0 =      A * ((A+1) - (A-1)*np.cos(w0) + 2 * np.sqrt(A) * alpha)
    b1 =  2 * A * ((A-1) - (A+1)*np.cos(w0))
    b2 =      A * ((A+1) - (A-1)*np.cos(w0) - 2 * np.sqrt(A) * alpha)
    a0 =          ((A+1) + (A-1)*np.cos(w0) + 2 * np.sqrt(A) * alpha)
    a1 = -2 *     ((A-1) + (A+1)*np.cos(w0))
    a2 =          ((A+1) + (A-1)*np.cos(w0) - 2 * np.sqrt(A) * alpha)

    b0 /=  a0
    b1 /=  a0
    b2 /=  a0
    a1 /=  a0
    a2 /=  a0
    
    coefs = np.array([np.array([b0, b1, b2]), np.array([a0, a1, a2])]);    
    return(coefs)
    
def sedea_rbj_hs(fc, fs, Q, gain):
    
    A = np.sqrt(10**(gain/20))
    w0 = 2 * np.pi * fc / fs
    alpha = np.sin(w0) / (2 * Q)
    
    b0 =      A *  ((A+1) - (A-1)*np.cos(w0) + 2 * np.sqrt(A) * alpha)
    b1 = -2 * A *  ((A-1) - (A+1)*np.cos(w0))
    b2 =      A *  ((A+1) - (A-1)*np.cos(w0) - 2 * np.sqrt(A) * alpha)
    a0 =           ((A+1) + (A-1)*np.cos(w0) + 2 * np.sqrt(A) * alpha)
    a1 =  2     *  ((A-1) + (A+1)*np.cos(w0))
    a2 =           ((A+1) + (A-1)*np.cos(w0) - 2 * np.sqrt(A) * alpha)

    b0 /=  a0
    b1 /=  a0
    b2 /=  a0
    a1 /=  a0
    a2 /=  a0
    
    coefs = np.array([np.array([b0, b1, b2]), np.array([a0, a1, a2])]);    
    return(coefs)

somenums = sedea_rbj_lpf(10000, 48000, 2.0, -1.0)
somenums1 = sedea_rbj_hpf(10000, 48000, 3.0, -1.0)
somenums2 = sedea_rbj_bpfcq(10000, 48000, 3.0, -1.0)
somenums3 = sedea_rbj_bpfcg(5000, 48000, 3.0, -1.0)
somenums4 = sedea_rbj_notch(5000, 48000, 3.0, -1.0)
somenums5 = sedea_rbj_apf(5000, 48000, 10, 1.0)
somenums6 = sedea_rbj_pek(5000, 48000, 10, 1.0)   
somenums7 = sedea_rbj_ls(5000, 48000, 1.0, -1.0)
somenums8 = sedea_rbj_hs(5000, 48000, 1.0, -1.0)          
#b = signal.firwin(80, 0.5, window=('kaiser', 8))

filterAnal(somenums[0], somenums[1], 48000)
filterAnal(somenums1[0], somenums1[1], 48000)
filterAnal(somenums2[0], somenums2[1], 48000)#
filterAnal(somenums3[0], somenums3[1], 48000)
filterAnal(somenums4[0], somenums4[1], 48000)
filterAnal(somenums5[0], somenums5[1], 48000)
filterAnal(somenums6[0], somenums6[1], 48000)
filterAnal(somenums7[0], somenums7[1], 48000)
filterAnal(somenums8[0], somenums8[1], 48000)
'''
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
'''