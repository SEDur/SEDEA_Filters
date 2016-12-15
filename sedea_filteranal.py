# -*- coding: utf-8 -*-
"""
Created on Fri Dec 09 09:57:05 2016

@author: simon
"""

def filterAnal(bcoefs, acoefs, SR):

    #b = signal.firwin(80, 0.5, window=('kaiser', 8))
    #b = np.array([somenums[0], somenums[1], somenums[2]])
    #a = np.array([somenums[3], somenums[4], somenums[5]])
    w, h = signal.freqz(bcoefs, acoefs)
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