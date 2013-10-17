
import numpy as np
import matplotlib.pyplot as plt
import random as rd

""" computation and visualisation of a discrete absorption-band model

    For a given array of target Q values, the code determines the optimal relaxation
    times and weights. This is done within in specified frequency range.


    C_r=1 and rho=1 are assumed. tau is computed from the target Q via 1/Q=0.5*tau*pi.
"""

#--------------------------------------------------------------------------------------------------
#- input ------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------

#- array of target Q's
Q_0=np.array([25.0, 100.0, 500.0, 1000.0])

#- minimum and maximum frequencies for optimisation in Hz
f_min=1.0/100.0
f_max=1.0/10.0

#- minimum and maximum frequencies for plotting in Hz
f_min_plot=1.0/200.0
f_max_plot=1.0/5.0

#- number of relaxation mechanisms
N=3

#- optimisation parameters (number of iterations, temperature, temperature decrease)
max_it=30000
T=0.2
d=0.9998

#--------------------------------------------------------------------------------------------------
#- initialisations --------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------

#- make logarithmic frequency axis
f=np.logspace(np.log10(f_min),np.log10(f_max),100)
w=2.0*np.pi*f

f_plot=np.logspace(np.log10(f_min_plot),np.log10(f_max_plot),100)
w_plot=2.0*np.pi*f_plot

#- compute tau from target Q
tau=2.0/(np.pi*Q_0)

#- compute initial relaxation times: logarithmically distributed relaxation times
tau_min=1.0/f_max
tau_max=1.0/f_min
tau_s=np.logspace(np.log10(tau_min),np.log10(tau_max),N)/(2*np.pi)

#- compute initial weights
D=1.0/tau_s
D=D/sum(D)

#- compute initial Q
chi=0.0

for n in np.arange(len(Q_0)):
    A=1.0
    B=0.0
    for p in np.arange(N):
        A+=tau[n]*(D[p]*w**2*tau_s[p]**2)/(1.0+w**2*tau_s[p]**2)
        B+=tau[n]*(D[p]*w*tau_s[p])/(1.0+w**2*tau_s[p]**2)

    Q=A/B
    chi+=sum((Q-Q_0[n])**2/Q_0[n]**2)

#--------------------------------------------------------------------------------------------------
#- search for optimal parameters ------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------

#- random search for optimal parameters
tau_s_test=np.array(np.arange(N),dtype=float)
D_test=np.array(np.arange(N),dtype=float)

for it in np.arange(max_it):

    #- compute perturbed parameters
    for p in np.arange(N):
        tau_s_test[p]=tau_s[p]*(1.0+(0.5-rd.random())*T)
        D_test[p]=D[p]*(1.0+(0.5-rd.random())*T)

    #- compute test Q
    chi_test=0.0

    for n in np.arange(len(Q_0)):
        A=1.0
        B=0.0
        for p in np.arange(N):
            A+=tau[n]*(D_test[p]*w**2*tau_s_test[p]**2)/(1.0+w**2*tau_s_test[p]**2)
            B+=tau[n]*(D_test[p]*w*tau_s_test[p])/(1.0+w**2*tau_s_test[p]**2)

        Q_test=A/B
        chi_test+=sum((Q_test-Q_0[n])**2/Q_0[n]**2)

    #- compute new temperature
    T=T*d

    #- check if the tested parameters are better
    if chi_test<chi:
        D[:]=D_test[:]
        tau_s[:]=tau_s_test[:]
        chi=chi_test

#--------------------------------------------------------------------------------------------------
#- print weights and relaxation times -------------------------------------------------------------
#--------------------------------------------------------------------------------------------------

print 'weights:  ', D
print 'relaxation times:  ', tau_s
print 'cumulative rms error:  ', np.sqrt(chi/(len(Q)*len(Q_0)))

#--------------------------------------------------------------------------------------------------
#- plot Q and phase velocity as function of frequency ---------------------------------------------
#--------------------------------------------------------------------------------------------------

for n in np.arange(len(Q_0)):

    #- compute optimal Q model for misfit calculations
    A=1.0
    B=0.0

    for p in np.arange(N):
        A+=tau[n]*(D[p]*w**2*tau_s[p]**2)/(1.0+w**2*tau_s[p]**2)
        B+=tau[n]*(D[p]*w*tau_s[p])/(1.0+w**2*tau_s[p]**2)

    Q=A/B
    chi=np.sqrt(sum((Q-Q_0[n])**2/Q_0[n]**2)/len(Q))

    print 'individual rms error for Q_0='+str(Q_0[n])+':  '+str(chi)

    #- compute optimal Q model for plotting
    A=1.0
    B=0.0

    for p in np.arange(N):
        A+=tau[n]*(D[p]*w_plot**2*tau_s[p]**2)/(1.0+w_plot**2*tau_s[p]**2)
        B+=tau[n]*(D[p]*w_plot*tau_s[p])/(1.0+w_plot**2*tau_s[p]**2)

    Q_plot=A/B
    v_plot=np.sqrt(2*(A**2+B**2)/(A+np.sqrt(A**2+B**2)))


    plt.subplot(121)
    plt.semilogx([f_min,f_min],[0.9*Q_0[n],1.1*Q_0[n]],'r')
    plt.semilogx([f_max,f_max],[0.9*Q_0[n],1.1*Q_0[n]],'r')
    plt.semilogx([f_min_plot,f_max_plot],[1.0*Q_0[n],1.0*Q_0[n]],'b')
    plt.semilogx(f_plot,Q_plot,'k',linewidth=2)
    plt.xlabel('frequency [Hz]')
    plt.ylabel('Q')
    plt.title('quality factor Q')

    plt.subplot(122)
    plt.semilogx([f_min,f_min],[0.9,1.1],'r')
    plt.semilogx([f_max,f_max],[0.9,1.1],'r')
    plt.semilogx(f_plot,v_plot,'k',linewidth=2)
    plt.xlabel('frequency [Hz]')
    plt.ylabel('v')
    plt.title('phase velocity')

    plt.show()

