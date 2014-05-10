# -*- coding: utf-8 -*-
import pylab as pl
import math
from scipy import integrate
import scipy as sp
import scipy.optimize
import numpy as np
h = [0, 0.25, 0.5, 1.0, 2.0, 4]
x = [0, 15, 30, 60, 120, 240]
y = [10, 9.178, 7.388, 6.247, 3.495, 1.38]


def trapezR(x, y):
    assert len(x) == len(y), 'input lists have to be of same length'
    trapezV = []
    for i in range(len(x)-1):
        currentValue = float((y[i]+y[i+1])/2)*(x[i+1]-x[i])
        trapezV.append(currentValue/60)
    return trapezV


def createReport(x,y, dose):
    estmatedAUC =  trapezR(x,y)
    hours = [float(i)/60 for i in x]
    for i in range(len(hours)-1):
        print 't = ' + str(hours[i]) + ': estimated AUC = ' + str(round(estmatedAUC[i], 2)) + ' mg*h/L'

    print 'Summed estimated AUC over timeframe = ' + str(round(sum(estmatedAUC), 2)) + ' mg*h/L'
    print 'Vd = ' + str(float(dose)/y[0])
    print 't(1/2) = ' + str(0.693*((float(dose)/y[0])/(float(dose)/sum(estmatedAUC))))
    print 'Clearance(estimated) = ' + str(round(float(dose)/sum(estmatedAUC),2))



def myExp(cStart, timeSteps, k):
    """
    takes a starting concentration and list of timesteps
    """
    estimate = []
    for time in range(timeSteps):
        estimate.append(cStart*np.exp(k*time))
    return estimate

def masterPlot(timeSteps, concentration):
    ### First make a List which contains all logged values of c
    logC = []
    for c in concentration:
        logC.append(np.log(c))
    ### compute k
    k, e = np.polyfit(x, logC, 1)
    ### get the estimates
    startC = concentration[0]
    endTime = timeSteps[-1]
    estTime = [i for i in range(endTime)]
    estConc = myExp(startC, endTime, k)
    ###
    exp = lambda x: startC*np.exp(k*x)
    AUC, err = integrate.quad(exp, 0, np.inf)
    text = 'AUC = ' + str(round(AUC/60))
    text2 = '$y= %d \cdot e^{-%rt}$' % (startC, round((k*60), 1))


    fig = pl.figure()
    fig.subplots_adjust(bottom=0.025, left=0.025, top = 0.975, right=0.975)

    pl.subplot(2, 1, 1)
    pl.scatter(timeSteps, concentration)
    pl.plot(estTime, estConc, color="green", label=text2)
    pl.title('Medikament X: AUC und approx. Formel beziehen sich auf Stunden')
    pl.ylabel('Konzentration c [mg/L]')
    pl.xlim(0, endTime)
    pl.grid()
    pl.legend()
    #pl.xticks(())
    #pl.yticks(())

    pl.subplot(2, 2, 3)
    pl.scatter(timeSteps, concentration)
    pl.ylabel('Linearisiert (log[c])')
    pl.grid()
    pl.xlim(0, endTime)
    pl.semilogy()
    pl.xlim(0)
    #pl.xticks(())
    #pl.yticks(())

    pl.subplot(2, 2, 4)
    pl.plot(estTime, estConc, label=text)
    pl.xlim(0, endTime)
    pl.ylabel('Regressionskurve')
    pl.fill_between(estTime, estConc, color='green', alpha=.25)
    pl.legend()
    #pl.xticks(())
    pl.yticks(())
    pl.grid()

    #pl.subplot(2, 3, 6)
    #pl.xticks(())
    #pl.yticks(())
    pl.show()

#masterPlot(x, y)
createReport(x,y,100)

def fit_exp_linear(t, y, C=0):
    y = y - C
    y = np.log(y)
    K, A_log = np.polyfit(t, y, 1)
    A = np.exp(A_log)
    return A, K

def model_func(t, A, K):
    return A * np.exp(K * t)

def fit_exp_nonlinear(t, y):
    opt_parms, parm_cov = sp.optimize.curve_fit(model_func, x, y, maxfev=1000)
    print opt_parms
    A, K = opt_parms
    return A, K
