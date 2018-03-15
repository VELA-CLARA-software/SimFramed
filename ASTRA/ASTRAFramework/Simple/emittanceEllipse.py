import numpy as np
from scipy.optimize import minimize
import math

def ellipse(x, y, a, b, alpha, x0, y0):
    return (((x - x0)*np.cos(alpha) + (y - y0)*np.sin(alpha))/a)**2 + (((x - x0)*np.sin(alpha) - (y - y0)*np.cos(alpha))/b)**2

def ellipseParticles(x, y, a, b, alpha, x0, y0):
    return zip(*[x,y,ellipse(x, y, a, b, alpha, x0, y0)])

def ellipseCut(x, y, a, b, alpha, x0, y0):
    allparticles = ellipse(x, y, a, b, alpha, x0, y0)
    return [i for i in allparticles if i < 1]

def ellipseCutParticles(x, y, a, b, alpha, x0, y0):
    allparticles = zip(*[x, y, ellipse(x, y, a, b, alpha, x0, y0)])
    return [i for i in allparticles if i[2] < 1]

def minFunc(inputvalues, x, y, confidence, normaliser):
    a, b, alpha = inputvalues
    survival = len(ellipseCut(x, y, abs(a), abs(b), abs(alpha), 0, 0))
    if survival < 1:
        survival = 1
    return 0.01*abs((math.pi*a*b)/normaliser)+(100*((float(survival)/float(len(x)))-confidence))**2

def emittanceEllipseOptimise(a, b, alpha, x, y, tolerance, confidence, iterations=5):
    x = x - np.mean(x)
    y = y - np.mean(y)
    simplexans = [(max(x) - min(x)), (max(y) - min(y)), 0]
    normaliser = confidence * math.pi * (max(x) - min(x)) * (max(y) - min(y))
    res = minimize(minFunc, simplexans, args=(x, y, confidence, normaliser),
            method='nelder-mead', options={'xtol': 1e-8, 'disp': False})
    a, b, c = res.x
    xans = 1.0/np.sqrt((np.cos(c)**2)/(a**2) + (np.sin(c)**2)/(b**2))
    yans = 1.0/np.sqrt((np.cos(c)**2)/(b**2) + (np.sin(c)**2)/(a**2))
    emit = abs(a*b)
    beta = emit / (yans**2)
    gamma = emit / (xans**2)
    alpha = np.sign(c) * np.sqrt((beta * gamma) - 1)
    return emit, alpha, beta, gamma
