''' Demo script to estimate pi using Monte Carlo sampling
'''
from numpy import random

def random_xy(xmin, xmax, ymin, ymax):
    x = random.uniform(xmin, xmax)
    y = random.uniform(ymin, ymax)
    return x, y

def check_circle(radius, x, y):
  return x**2 + y**2 < radius**2

radius = 1.
NTRIALS = 1000000000
random.seed()

in_circle = 0
for n in range(NTRIALS):
  new_x, new_y = random_xy(-radius, radius, -radius, radius)
  if check_circle(radius, new_x, new_y):
    in_circle += 1
  if n % 1000000 == 0:
    print(n, in_circle / (n + 1) * 4.)
