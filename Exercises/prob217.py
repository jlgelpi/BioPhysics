#! /usr/bin/python

__author__ = "gelpi"
__date__ = "$02-nov-2017 8:21:31$"

import math

#Parameters

SIG = 3.4
EPS = 0.09


class System():

    """Represents the ensemble of particles"""

    def __init__(self, dmin):
        """Generates a regular lattice of 3x3x3 particles in internal coordinates separated by dmin"""
        self.parts = []
        self.d = dmin
        for x in range(-1, 2):
            for y in range(-1, 2):
                for z in range(-1, 2):
                    self.parts.append(Particle(x, y, z, 0.))

    def calc_energy(self):
        """Calculates vdw + elec energy for the whole system"""
        eint = 0.
        evdw = 0.
        for p1 in self.parts:
            for p2 in self.parts:
                if p1 == p2:
                    continue
                eint = eint + 0.5 * p1.elec_interaction(p2, self.d)
                evdw = evdw + 0.5 * p1.vdw_interaction(p2, self.d)
        return [evdw, eint]


class Particle():
    """ A single particle,has coordinates (in internal units) and charge (in e) """
    def __init__(self, x, y, z, c):
        self.x = x
        self.y = y
        self.z = z
        self.c = c

    def distance(self, other):
        """Evaluates Sqrt[(x0-x1)2 +(y0-y1)2 + (z0-z1)2]"""
        return math.sqrt((self.x-other.x) ** 2 + (self.y-other.y) ** 2 + (self.z-other.z) ** 2)

    def vdw_interaction(self, other, dmin):
        """vdwenergy between two particles"""
        f = SIG / (self.distance(other) * dmin)
        return 4. * EPS * (pow(f, 12)-pow(f, 6))

    def elec_interaction(self, other, dmin):
        """Electrostatic interaction"""
        d = self.distance(other) * dmin
        return 332.16 * self.c * other.c / d

    def __eq__(self, other):
        """Allow to decide whether two particles are the same"""
        return self.x == other.x and self.y == other.y and self.z == other.z

sys = System(3.8) #dmin included in the code, should be an external parameter!!

[evdw, eint] = sys.calc_energy()

print ("Evdw=", evdw)
#Removing central particle
#p0 is the central particle
p0 = Particle(0, 0, 0, 0)
evdw0 = 0
for pi in sys.parts:
    if pi != p0:
        evdw0 = evdw0 + p0.vdw_interaction(pi, sys.d)
print ("New Evdw=", evdw-evdw0)

#Electrostatics, equal positive charge
#Calculate for a q=1e and adjust later
for pi in sys.parts:
    pi.c = 1.
(evdw1, eint1) = sys.calc_energy()
print ("Charge to equilibrate: ",math.sqrt(-evdw / eint1))

# Negative in central particle
ecen=0.
p0.c=1.
for pi in sys.parts:
    if pi != p0:
        ecen = ecen +  p0.elec_interaction(pi,sys.d)
neweint = eint1 - 2. * ecen
print ("Centrl neg to compensate: ", math.sqrt(-evdw/neweint))

