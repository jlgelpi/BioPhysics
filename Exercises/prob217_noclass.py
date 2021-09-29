#! /usr/bin/python

__author__ = "gelpi"
__date__ = "$02-nov-2017 8:21:31$"

import math

#Parameters

SIG = 3.4
EPS = 0.09
dmin = 3.8

def calc_energy(sys):
    """Calculates vdw + elec energy for the whole system"""
    eint = 0.
    evdw = 0.
    for p1 in sys['parts']:
        for p2 in sys['parts']:
            if p1 == p2:
                continue
            eint = eint + 0.5 * elec_interaction(p1, p2, sys['d'])
            evdw = evdw + 0.5 * vdw_interaction(p1, p2, sys['d'])
    return [evdw, eint]

def distance(p1,p2):
    """Evaluates Sqrt[(x0-x1)2 +(y0-y1)2 + (z0-z1)2]"""
    return math.sqrt((p1[0]-p2[0]) ** 2 + (p1[1]-p2[1]) ** 2 + (p1[2]-p2[2]) ** 2)

def vdw_interaction(p1, p2, dmin):
    """vdwenergy between two particles"""
    f = SIG / (distance(p1,p2) * dmin)
    return 4. * EPS * (pow(f, 12)-pow(f, 6))

def elec_interaction(p1,p2, dmin):
    """Electrostatic interaction"""
    d = distance(p1,p2) * dmin
    return 332.16 * p1[3] * p2[3] / d

# MAIN 

sys = {'parts': [], 'd' : dmin}

for x in range(-1, 2):
    for y in range(-1, 2):
        for z in range(-1, 2):
            sys['parts'].append([x, y, z, 0.])

[evdw, eint] = calc_energy(sys)

print ("Evdw=", evdw)
#Removing central particle
#p0 is the central particle
p0 = [0, 0, 0, 0]
evdw0 = 0
for pi in sys['parts']:
    if pi != p0:
        evdw0 = evdw0 + vdw_interaction(p0, pi, sys['d'])
print ("New Evdw=", evdw-evdw0)

#Electrostatics, equal positive charge
#Calculate for a q=1e and adjust later
for pi in sys['parts']:
    pi[3] = 1.
(evdw1, eint1) = calc_energy(sys)
print ("Charge to equilibrate: ",math.sqrt(-evdw / eint1))

# Negative in central particle
ecen=0.
p0[3]=1.
for pi in sys['parts']:
    if pi != p0:
        ecen = ecen +  elec_interaction(p0,pi,sys['d'])
neweint = eint1 - 2. * ecen
print ("Centrl neg to compensate: ", math.sqrt(-evdw/neweint))

