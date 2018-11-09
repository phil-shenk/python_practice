import numpy as np
import csv
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import scipy.signal
import os
import operator

def line(x, m, b):
    return m*x+b

def diode_curve(V, I0, V0, b, m):
    return I0*(np.exp(V/V0) - 1)+b+(V*m)
    
def getFlatFit(Varr, Iarr, points):
    ####fit a line to the first segment
    flatFitPoints = points  #number of points to use for flat part of line fit (out of total 200)
    #do the line fit on current vs voltage
    fit_params, fit_uncertainties = curve_fit(line, Varr[0:flatFitPoints], Iarr[0:flatFitPoints])
    slope = fit_params[0]
    intercept = fit_params[1]
    #tangent line to flat section
    flatFit = slope*np.array(Varr) + intercept
    return flatFit
    
def iterativeDerivative(Varr, Iarr):
    dIdV = []
    for i in range(len(Varr)-1):
        slope = (Iarr[i+1] - Iarr[i]) / (Varr[i+1] - Varr[i])
        dIdV.append(slope) 
    return dIdV
        
        
def findMax(Xarr, Yarr):
    # find max slope
    maxIndex, maxValue = max(enumerate(Yarr), key=operator.itemgetter(1))
    return maxIndex, maxValue
    
def thresholdProxy(Varr, smoothI, flatFit, threshold):
    #### STOPPING POTENTIAL PROXY METHOD 1: Threshold
    #threshold check, go along and compare height to flat fit line
    thresholdHit = False
    exceedIndex = 0
    for i in range(len(Varr)):
        difference = smoothI[i] - flatFit[i]
        if difference > threshold and not thresholdHit:
            thresholdHit = True
            exceedIndex = i    
    #exceedIndex = where the (smoothed) current rises above flatFit+threshold
    return exceedIndex
    
def doubleIntersectionProxy(Varr, slopeFit, flatFit):
    #### STOPPING POTENTIAL PROXY METHOD 2: Tangent Intsersection
    #tangent intersection method, find where the flat fit and steepest slope fit intercect   
    intersectionHit = False
    intersectionIndex = 0
    for i in range(len(Varr)):
        difference = slopeFit[i] - flatFit[i]
        if(difference > 0) and not intersectionHit:
            intersectionHit = True
            intersectionIndex = i
    #intersectionIndex = where the flat tangent and sloped tangent intersect
    return intersectionIndex

def curveFitProxy(Varr, Iarr, maxSlopeIndex, flatFit, i):
    #### STOPPING POTENTIAL PROXY METHOD 3: Fitting a diode curve
    #do the line fit on current vs voltage
    I0 = 0
    V0 = 1
    b = 0
    m = 0
    
    #subtract flatfit from I data to align along 0
    #correctedIarr = Iarr #- flatFit    
    
    if(maxSlopeIndex > 0):
        fit_params, fit_uncertainties = curve_fit(diode_curve, Varr[0:maxSlopeIndex], Iarr[0:maxSlopeIndex], maxfev = 50000)
        I0 = fit_params[0]
        V0 = fit_params[1]
        b  = fit_params[2]
        m  = fit_params[3]

    print('Fitted V0 value for '+str(i)+' is '+str(V0))
    diodeFit = I0*(np.exp(np.array(Varr)/V0) - 1)+b+(m*np.array(Varr))
    
    plt.figure(100+fileIndex)
    plt.ylim(-1,5)
    plt.plot(Varr[0:len(Varr)-0], diodeFit)
    plt.scatter(Varr, Iarr)
    
    return diodeFit, V0 
    
def xinterceptProxy(Varr, smoothI):
    #### STOPPING POTENTIAL PROXY METHOD 4: X intercept
    #find where the smoothed current crosses the x intercept, nice and simple
    intersectionHit = False
    intersectionIndex = 0
    for i in range(len(smoothI)):
        if(smoothI[i] > 0) and not intersectionHit:
            intersectionHit = True
            intersectionIndex = i
    #intersectionIndex = where the flat tangent and sloped tangent intersect
    return intersectionIndex
    
def AnalyzeCurve(Varr, Iarr, number, intensity, previousFigure):
    
    # Calculate flat fit line for first flat section (first 50 points works)
    flatFit = getFlatFit(Varr, Iarr, 50)    
    
    # Smooth current with Savitsky-Golay filter
    smoothI = scipy.signal.savgol_filter(Iarr, 11, 3)
    
    # take derivative of smoothed current
    smoothdIdV = iterativeDerivative(Varr, smoothI)
    
    # find max slope
    maxSlopeIndex, maxSlopeValue = findMax(Varr, smoothdIdV)
    
    #find ordered pair of the max slope point (for point slope form)
    xpos = Varr[maxSlopeIndex]
    ypos = smoothI[maxSlopeIndex]
    
    #y-y1 = m(x-x1) => y = m(x-x1)+y1 where x1,y1 is a point on the IvsV graph
    #tangent line to slope section
    slopeFit = maxSlopeValue*(np.array(Varr)-xpos)+ypos
    
    #### STOPPING POTENTIAL PROXY METHOD 1: Threshold
    exceedIndex = thresholdProxy(Varr, smoothI, flatFit, 0.025)
    
    #### STOPPING POTENTIAL PROXY METHOD 2: Tangent Intsersection
    intersectionIndex = doubleIntersectionProxy(Varr, slopeFit, flatFit)
            
    #### STOPPING POTENTIAL PROXY METHOD 3: Curve Fit            
    diodeFit, V0 = curveFitProxy(Varr, Iarr, maxSlopeIndex, flatFit, number)
    plt.figure(previousFigure)
    
    #### STOPPING POTENTIAL PROXY METHOD 4: X intercept
    xinterceptV0 = xinterceptProxy(Varr, smoothI)
            
    #return both fit lines, smoothed Current and Derivative of Current, and both proxy indices
    return flatFit, smoothI, smoothdIdV, slopeFit, diodeFit, exceedIndex, intersectionIndex, V0, xinterceptV0
    
#which file is being read
fileIndex=0
#values of stopping potential obtained with each proxy method
stoppingPotentials_threshold=[] #calculated with threshold method
stoppingPotentials_intercept=[] #calculated with intercept method
stoppingPotentials_diodefit=[] #calculated with intercept method
stoppingPotentials_xintFit=[]
frequencies=[]
colorFrequencies = {
'blue.dat':435.8,
'green.dat':546.1,
'yellow.dat':577.0,
'violet.dat':404.7,
'ultraviolet.dat':365.4
}
for filename in os.listdir('data'):
    rows = []
    print('data/'+filename)
    with open('data/'+filename) as f:
        rows = f.readlines()
    #print(rows,'rows')
        
    dataRows = []
    voltageArrays=[]
    currentArrays=[]
    voltBuffer = []
    currentBuffer = []
    
    for textline in rows:    #skip lines that aren't data
        #print(line)
        if (len(textline.strip()) == 0 or textline.startswith('Time') or textline.startswith('Voltage')):
            continue
        #dump array into big array if you are done with the array
        elif(textline.startswith('Run')):
            voltageArrays.append(list(voltBuffer))
            voltBuffer = []
            currentArrays.append(list(currentBuffer))
            currentBuffer = []
        else:
            #print('case3 ***' +line)
            vals = textline.split()
            floatvals = []
            voltBuffer.append(float(vals[0]))
            currentBuffer.append(float(vals[1]))
    voltageArrays.append(list(voltBuffer))
    voltBuffer = []
    currentArrays.append(list(currentBuffer))
    currentBuffer = []
        
    currentArrays.remove(currentArrays[0])
    voltageArrays.remove(voltageArrays[0])
    #print(voltageArrays)
    #print(currentArrays)
    
    
    #plot it
    plt.title("Data from " +filename)
    plt.ylabel("Y")
    plt.xlabel("X")
    
    
    percentages = [0,20,40,60,80,100,0]
    for i in range(len(voltageArrays)):
        plt.figure(fileIndex)
        plt.xlabel('V')
        plt.ylabel('I')
        plt.ylim(-1,5)
        plt.title('V vs I with fits, '+filename)
        #plt.plot(voltageArrays[i], currentArrays[i], label=(str(i)+': '+str(percentages)+'%'))
        flatFit, smoothI, smoothdIdV, slopeFit, diodeFit, exceedIndex, intersectionIndex, diodeFitV0, xintV0 = AnalyzeCurve(voltageArrays[i], currentArrays[i], i, percentages[i], fileIndex)
        plt.subplot(211)        
        plt.plot(voltageArrays[i], flatFit)
        plt.plot(voltageArrays[i], smoothI, label=str(percentages[i])+'%')
        plt.subplot(212)
        plt.plot(voltageArrays[i], flatFit)
        plt.plot(voltageArrays[i], smoothI, label=str(percentages[i])+'%')
        
        # log stopping potentials from threshold method
        plt.subplot(212)
        plt.legend(loc='best')
        plt.scatter(voltageArrays[i][exceedIndex], smoothI[exceedIndex])
        if not (percentages[i] == 0):
            stoppingPotential = voltageArrays[i][exceedIndex]
            stoppingPotentials_threshold.append(stoppingPotential)
            frequencies.append(colorFrequencies.get(filename))
        
        plt.subplot(211)
       
        #log stopping potentials from intersection method
        plt.scatter(voltageArrays[i][intersectionIndex], flatFit[intersectionIndex])
        if not (percentages[i] == 0):
            stoppingPotential = voltageArrays[i][intersectionIndex]
            stoppingPotentials_intercept.append(stoppingPotential)
            
            #also add diode fit V0
            stoppingPotentials_diodefit.append(diodeFitV0)
            #also xintercept
            stoppingPotentials_xintFit.append(xintV0)
        
        plt.subplot(211)
        plt.plot(voltageArrays[i], slopeFit)
        #plt.figure(1)
        #plt.xlabel('V')
        #plt.ylabel('smooth dI/dV')
        #plt.title('Derivative of Smoothed Current')
        #plt.plot(voltageArrays[i][0:len(voltageArrays[i])-1], smoothdIdV)
        plt.legend(loc='best')
   

    fileIndex = fileIndex+1
    
plt.figure(25)
plt.xlabel('wavelength')
plt.ylabel('stopping potential (V)')
plt.scatter(frequencies, stoppingPotentials_intercept, c='b', marker='*', s=100)
plt.scatter(frequencies, stoppingPotentials_threshold, c='r')
plt.scatter(frequencies, stoppingPotentials_diodefit, c='y')
plt.scatter(frequencies, stoppingPotentials_xintFit, c='g')
#plt.plot(voltageArrays,currentArrays)
plt.show()
