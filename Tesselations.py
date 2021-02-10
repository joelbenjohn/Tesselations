import itertools
import numpy
from scipy.spatial import ConvexHull
import math
from math import pi
from matplotlib.collections import LineCollection
from matplotlib import pyplot as plot
from copy import deepcopy
from functools import reduce
import pdb
from zipfile import ZipFile 
import os 
def rectangle(n1, n2, Side):
  S = numpy.zeros(((n1+2)*(n2+2), 2))
  R = numpy.ones((n1+2)*(n2+2))
  Side = 5; x = Side/(n1); y = Side/(n2);
  lim = numpy.array([[0,0], [Side, Side]])
  for i in range(n1+2):
      for j in range(n2+2):
              S[(n2+2)*i+j, 0] = lim[0, 0]- x/2 + x*j
              S[(n2+2)*i+j, 1] = lim[0, 1]- y/2 + y*i
  return S, R

def triangle(Seeds, Side, Theta = numpy.sqrt(3)):
  n1 = int(numpy.sqrt(3*Seeds))+1
  n2 = int(n1/2.8)+1
  Seeds = n1*n2
  S = numpy.zeros(((n1+4)*(n2+4), 2))
  R = numpy.ones((n1+4)*(n2+4))
  z = Side/(2*Theta*Seeds)*(-1+numpy.sqrt(1+4*Theta*Seeds))+Side/200
  y = z*Theta
  x = (Theta**2-1)/(2*Theta)*z
  for i in range(n1+4):
      for j in range(n2+4):
          if((i+1)%4==1):
              S[(n2+4)*i+j, 0] = 2*z*j - 4*z
              S[(n2+4)*i+j, 1] = x + 2*y*math.floor(i/4) - 4*y

          if((i+1)%4==2):
              S[(n2+4)*i+j, 0] = z + 2*z*j- 4*z
              S[(n2+4)*i+j, 1] = y-x + 2*y*math.floor(i/4) - 4*y

          if((i+1)%4==3):
              S[(n2+4)*i+j, 0] = z + 2*z*j- 4*z
              S[(n2+4)*i+j, 1] = x + y +2*y*math.floor(i/4) - 4*y

          if((i+1)%4==0):
              S[(n2+4)*i+j, 0] = 2*z*j- 4*z
              S[(n2+4)*i+j, 1] = y-x + y + 2*y*math.floor(i/4) - 4*y
  return S, R

def honeycomb(Seeds, Side):
  x = Side/numpy.sqrt(2*Seeds*numpy.sqrt(3))
  n1 = int(Side/x/numpy.sqrt(3))
  n2 = int(Seeds/n1)+1
  R = numpy.ones((n1+100)*(n2+100))
  S = numpy.zeros(((n1+100)*(n2+100), 2))
  for i in range(n1+100):
      for j in range(n2+100):
          if((j+1)%2==1):
              S[(n2+100)*i+j, 0] = numpy.sqrt(3)*x*j-25*numpy.sqrt(3)
              S[(n2+100)*i+j, 1] = 2*x*i - 50*x

          if((j+1)%2==0):
              S[(n2+100)*i+j, 0] = numpy.sqrt(3)*x*j - 25*numpy.sqrt(3)
              S[(n2+100)*i+j, 1] = x + 2*x*i - 50*x
  return S, R

def semi1(Seeds, Side, w):
  r = Side*numpy.sqrt(2/Seeds/numpy.sqrt(3))
  n2 = int(Side/r)+1
  n1 = int(Side*2/r/numpy.sqrt(3))+1
  R = numpy.ones((n1+8)*(n2+8))
  S = numpy.zeros(((n1+8)*(n2+8), 2))
  # x = 1
  # r = x*(0.5/numpy.tan(pi/12)+1/numpy.sqrt(3))
  wo = r*w/(1+w)
  wt = r/(1+w)
  for i in range(n2+8):
      for j in range(n1+8):
          if (j+1)%2 == 1:
              S[(n1+8)*i+j, 0] = r*numpy.sin(pi/3)*j - 4*r*numpy.sin(pi/3)
              S[(n1+8)*i+j, 1] = r*i - 4*r
              if(i+1)%3 == 1 or (i+1)%3 == 0:
                  R[(n1+8)*i+j] = wt
              if(i+1)%3==2:
                  R[(n1+8)*i+j] = wo
          if (j+1)%2 == 0:
              S[(n1+8)*i+j, 0] = r*numpy.sin(pi/3)*j -4*r*numpy.sin(pi/3)
              S[(n1+8)*i+j, 1] = 0.5*r + r*i - 4*r
              if(i+1)%3==1 or (i+1)%3==2:
                  R[(n1+8)*i+j] = wt
              if(i+1)%3==0:
                  R[(n1+8)*i+j] = wo
  return S, R

def pert(Seeds, Side, w, c, g):
  r = Side*numpy.sqrt(2/Seeds/numpy.sqrt(3))
  n2 = int(Side/r)+1
  n1 = int(Side*2/r/numpy.sqrt(3))+1
  R = numpy.empty((0))
  S = numpy.empty((0, 2))
  # x = 1
  # r = x*(0.5/numpy.tan(pi/12)+1/numpy.sqrt(3))
  wo = r*w/(1+w)/2
  wt = r/(1+w)/2
  for i in range(n2+8):
      for j in range(n1+8):
          if (j+1)%2 == 1:
            count = 0
            attempt = 0 
            while (count == 0):
              s0 = r*numpy.sin(pi/3)*j - 4*r*numpy.sin(pi/3) + numpy.random.normal(0, c*r/32)
              s1 = r*i - 4*r + numpy.random.normal(0, c*r/32)
              if (i+1)%3 == 1 or (i+1)%3 == 0:
                  r0 = numpy.abs(numpy.random.normal(wt, g*wt/64))
              if (i+1)%3 == 2:
                  r0 = numpy.abs(numpy.random.normal(wo, g*wo/64))
              if len(S) == 0:
                S = numpy.append(S, [[s0, s1]], axis = 0)
                R = numpy.append(R, [r0], axis = 0)
              elif numpy.all(numpy.sum((S - [[s0, s1]])**2, axis = 1)-(R+r0)**2>=-0.6*r):
                S = numpy.append(S, [[s0, s1]], axis = 0)
                R = numpy.append(R, [r0], axis = 0)   
                count += 1  
              if attempt>500:
                S = numpy.append(S, [[r*numpy.sin(pi/3)*j - 4*r*numpy.sin(pi/3), r*i - 4*r]], axis = 0)
                R = numpy.append(R, [wo/4] if (i+1)%3 == 2 else [wt/4], axis = 0)
                count += 1
              attempt += 1          
          if (j+1)%2 == 0:
            count = 0
            attempt = 0 
            while (count == 0):
              s0 = r*numpy.sin(pi/3)*j - 4*r*numpy.sin(pi/3) + numpy.random.normal(0, c*r/32)
              s1 = 0.5*r + r*i - 4*r + numpy.random.normal(0, c*r/32)
              if (i+1)%3 == 1 or (i+1)%3 == 2:
                  r0 = numpy.abs(numpy.random.normal(wt, g*wt/64))
              if (i+1)%3 == 0:
                  r0 = numpy.abs(numpy.random.normal(wo, g*wo/64))
              if len(S) == 0:
                S = numpy.append(S, [[s0, s1]], axis = 0)
                R = numpy.append(R, [r0], axis = 0)
              elif numpy.all(numpy.sum((S - [[s0, s1]])**2, axis = 1)-(R+r0)**2>=-0.6*r):
                S = numpy.append(S, [[s0, s1]], axis = 0)
                R = numpy.append(R, [r0], axis = 0)
                count += 1
              if attempt>500:
                S = numpy.append(S, [[r*numpy.sin(pi/3)*j - 4*r*numpy.sin(pi/3), 0.5*r + r*i - 4*r]], axis = 0)
                R = numpy.append(R, [wo/4] if (i+1)%3 == 2 else [wt/4], axis = 0)
                count += 1
              attempt += 1
  return S, R
def semi2(Seeds, Side, w):
  r = Side/numpy.sqrt(Seeds)
  n2 = int(numpy.sqrt(Seeds)+1)
  n1 = int(numpy.sqrt(Seeds)+1)
  R = numpy.ones((n1+8)*(n2+8))
  S = numpy.zeros(((n1+8)*(n2+8), 2))
  # x = 1
  # r = x*(0.5/numpy.tan(pi/12)+1/numpy.sqrt(3))
  wo = r*w/(1+w)
  wt = r/(1+w)
  for i in range(n2+8):
      for j in range(n1+8):
          S[(n1+8)*i+j, 0] = r*j - 4*r
          S[(n1+8)*i+j, 1] = r*i - 4*r
          if(i+j)%2 == 1:
              R[(n1+8)*i+j] = wt
          if(i+j)%2==0:
              R[(n1+8)*i+j] = wo
  return S, R

def semi3(Seeds, Side, w):
  r = Side*numpy.sqrt(2/Seeds/numpy.sqrt(3))/1.5
  n2 = int(Side/r)+1
  n1 = int(Side*2/r/numpy.sqrt(3))+1
  f = numpy.array([1, 1, 2, 2])
  f2 = numpy.array([0, 0, 1, 0])
  f1 = numpy.zeros(n1+16)
  for i in range(1, n1+16):
    f1[i] = f1[i-1]+ f[(i-1)%4]
  ws = r
  r3 = numpy.sqrt(3)
  wt = r*(numpy.sqrt(1-2/3/(1/r3+1)**2))
  wh = r*(numpy.sqrt(1+2/(1/r3+1)**2))
  S = []
  R = []
  for i in range(n2+16):
      if(i+1)%2 == 1:
        for j in range(n1+16):
              S.append([f1[j]*r-4*r+ f2[(i)%4]*(2*r+2*r*numpy.cos(pi/3)), r*i*numpy.sin(pi/3) - 4*r*numpy.sin(pi/3)])
              if(j+1)%4 == 1 or (j+1)%4 == 3:
                  R.append(wt)
              if(j+1)%4==2:
                  R.append(ws)
              if(j+1)%4==0:
                  R.append(wh)
      if(i+1)%2 == 0:
          for j in range(int((n1+16)/2)):
              S.append([2*r*(numpy.cos(pi/3)+1)*j -4*r-r*numpy.cos(pi/3), r*i*numpy.sin(pi/3) - 4*r*numpy.sin(pi/3)])
              R.append(ws)
  S = numpy.array(S)
  R = numpy.array(R)
  return S, R

def semi4(Seeds, Side, w):
  r = Side*numpy.sqrt(2/Seeds/numpy.sqrt(3))
  n2 = int(Side/r)+1
  n1 = int(Side*2/r/numpy.sqrt(3))+1
  f = [0,-1]
  S=[];R=[]
  t= 1/numpy.tan(pi/6)
  t1 = 1/numpy.tan(pi/12)
  ws =r; wh = r*numpy.sqrt(((t-1)/3/(t+1)+1)); wo = r*numpy.sqrt(((t1-1)/(t1+1)+1))
  f1 = [ws, wo]
  f2 = [r/numpy.sqrt(3), r*numpy.sqrt(4/3)]
  for i in range(4):
    if (i+1)%4==1:
      for j in range(n1+8):
        S.append([-4*r+r*j+f[0 if i%8==0 else 1]*r, -numpy.sqrt(3)+numpy.sqrt(3)*r*i/4])
        R.append(f1[j%2])
    if (i+1)%4==2 or (i+1)%4==0:
      for j in range(int((n1+8)/2)):
        S.append([-4*r+2*r*j+f[0 if (i-1)%4==0 else 1]*r+f[0 if ((i-1)%8==0 or (i+1)%8==0) else 1]*2*r, -numpy.sqrt(3)+numpy.sqrt(3)*r*numpy.floor(i/4)+f2[0 if (i-1)%4==0 else 1]])
        R.append(wh)
    if (i+1)%4==3:
      for j in range(n1+8):
        S.append([-4*r+r/2+r*j+f[0 if (i-2)%8==0 else 1]*r, -numpy.sqrt(3)+r*numpy.sqrt(3)/2+numpy.sqrt(3)*r*numpy.floor(i/4)])
        R.append(ws)
  S = numpy.array(S)
  R = numpy.array(R)
  S0 = deepcopy(S)
  R0 = deepcopy(R)
  for i in range(int(n2/1.5)):
    R = numpy.append(R, R0)
    if (i%2==0):
      S1 = S0+numpy.array([-r, numpy.sqrt(3)*r*(i+1)])
    else:
      S1 = S0+numpy.array([0, numpy.sqrt(3)*r*(i+1)])
    S = numpy.append(S, S1, axis=0)
    S = numpy.array(S)
  R = numpy.array(R)
  return S, R

def semi5(Seeds, Side, w):
  r = Side*numpy.sqrt(2/Seeds/numpy.sqrt(3))
  n2 = int(Side/r)+1
  n1 = int(Side*2/r/numpy.sqrt(3))+1
  R = numpy.ones((n1+8)*(n2+8))
  S = numpy.zeros(((n1+8)*(n2+8), 2))
  # x = 1
  # r = x*(0.5/numpy.tan(pi/12)+1/numpy.sqrt(3))
  t = 1/numpy.tan(pi/8)
  wo = r*numpy.sqrt(((t-1)/(t+1)+1))
  wt = r
  for i in range(n2+8):
      for j in range(n1+8):
          S[(n1+8)*i+j, 0] = r*j - 4*r
          S[(n1+8)*i+j, 1] = r*i - 4*r
          if(i+j)%2 == 1:
              R[(n1+8)*i+j] = wt
          if(i+j)%2==0:
              R[(n1+8)*i+j] = wo
  return S, R
  
def semi6(Seeds, Side, w):
  r = Side*numpy.sqrt(2/Seeds/numpy.sqrt(3))
  n2 = int(Side/r)+1
  n1 = int(Side*2/r/numpy.sqrt(3))+1
  t = 1/numpy.tan(pi/3)
  wt = r*numpy.sqrt(((t-1)/(t+1)+1))
  ws = r
  S=[];R=[]
  for i in range(3):
      for j in range(n1+8):
          S.append([2*r/(1+t)*j - 4*r, r*i - 4*r])
          if (i)%3 == 0 or i%3==2:
              R.append(wt)
          if(i)%3==1:
              R.append(ws)
  S = numpy.array(S)
  R = numpy.array(R)
  S0 = deepcopy(S)
  R0 = deepcopy(R)
  for i in range(int(n2/2)):
    R = numpy.append(R, R0)
    if (i%2==0):
      S1 = S0+numpy.array([-r/(1+t), r*(2*t+1)/(t+1)/t*(i+1)])
    else:
      S1 = S0+numpy.array([0, r*(2*t+1)/(t+1)/t*(i+1)])
    S = numpy.append(S, S1, axis=0)
    S = numpy.array(S)
  R = numpy.array(R)
  return S, R

def mono(Seeds, Side, k):
  R = numpy.ones(int(Seeds*1.2))
  # S = numpy.random.uniform(-Side/5, 6*Side/5, (int(Seeds*1.2), 2))
  S = numpy.empty((0, 2))
  # S = numpy.array([[Side/5, Side/5], [Side/5, Side/5+6], [Side/5+6, Side/5], [Side/5+6, Side/5+6]])
  R = numpy.empty((0))
  # R = numpy.array([3, 3, 3, 3])
  Smax = [-Side/5, -Side/5.01]
  Smean =  [0, 0.1]
  count = 0
  m = numpy.sqrt(Side**2/Seeds/numpy.pi)
  R1 = numpy.abs(numpy.array(numpy.random.normal(m, k, (int(1.96*Seeds*1.5)))))
  R11 = numpy.sort(R1)[::-1]
  count1= 0
  R0 = [R11[0]]
  nseed = 0
  while len(S) < int(Seeds*1.96):
    count1+=1
    # Smax = numpy.max(S, axis = 0)
    # Smean = numpy.mean(S, axis = 0)
    # R0 = numpy.random.rand(1)*7
    # R0 = numpy.array([3])

    if len(S)==0:
      S0 = [[-Side/5, -Side/5.01]]
      S = numpy.append(S, S0, axis = 0)
      R= numpy.append(R, R0, axis = 0)
    # if Smean[0]<Smean[1]:
    #   S01 = numpy.random.uniform(-Side/5, Smax[1])
    #   # print(S - [2*Smean[0], S01])
    #   # print(numpy.sum((S - [2*Smean[0], S01])**2, axis = 1))
    #   ind = numpy.argmin(numpy.sum((S - [Side, S01])**2, axis = 1))
    #   # print(S[ind, 0], R0[0], R[ind])
    #   S00 = S[ind, 0]+ 2*R0[0] + R[ind]
    #   # print(S00)
    # else:
    #   S00 = numpy.random.uniform(-Side/5, Smax[0])
    #   ind = numpy.argmin(numpy.sum((S - [Side, S00])**2, axis = 1))
    #   # print(S[ind, 1], R0[0], R[ind])
    #   S01 = S[ind, 1]+ 2*R0[0] + R[ind]
    #   # print(S01)
    # S00 = -Side/5 + (1.8+(k)/25)*Smean[0] + numpy.random.uniform(0, 2*R0[0]) if Smean[0]<Smean[1] else numpy.random.uniform(-Side/5, Smax[0])
    # S01 = -Side/5.01 + (1.8+k/25)*Smean[1] + numpy.random.uniform(0, 2*R0[0]) if Smean[1]<=Smean[0] else numpy.random.uniform(-Side/5, Smax[1])
    # S0 = [[S00, S01]]
    S0 = numpy.random.uniform(-Side/5, Side*8/5, (1, 2))
    if count1>200:
      # R0 = numpy.abs(numpy.array([numpy.random.normal(5, 2)]))
      if nseed+1<int(1.96*Seeds*1.5):
        nseed += 1
        R0 = [R11[nseed]]
      else:
        Seeds = Seeds*0.9
    if numpy.all(numpy.sum((S-S0)**2, axis = 1)-R**2-R0[0]**2>0):
      Smean = numpy.mean(S, axis = 0) - [-Side/5, -Side/5.01] + 2
      Smax = numpy.max(S, axis = 0)
      S = numpy.append(S, S0, axis = 0)
      R= numpy.append(R, R0, axis = 0)
      # R0 = numpy.abs(numpy.array([numpy.random.normal(5, 2)]))
      if nseed+1<int(1.96*Seeds*1.5):
        nseed += 1
        R0 = [R11[nseed]]
      else:
        Seeds = Seeds*0.9
      count1 = 0
    count+=1
    
  return S, R

def poly(Seeds, Side, Centers, Size):
  R = numpy.random.beta(1, 1000, int(Seeds*1.2))*5
  S1 = numpy.random.rand(Centers, 2)*Side
  S0 = (numpy.random.rand(int(Seeds*1.2), 2)*1.4-0.2)*Side
  d=[[] for i in range(Centers)]
  for i in range(len(S1)):
    di = numpy.sum((S0-S1[i])**2, axis = 1)-Size**2
    index = numpy.argwhere(di>0)
    d[i] = index
  index = reduce(numpy.intersect1d, (d[i] for i in range(Centers)))
  S = S0[index, :]
  R = numpy.ones(len(index))
  siz = (numpy.random.rand(Centers)+1)*Size
  S=numpy.append(S, S1, axis = 0)
  R=numpy.append(R, siz, axis = 0)
  return S, R

def voronoi1(Seeds, Side, dev):
  x = Side/numpy.sqrt(Seeds)
  S = numpy.empty((0,2))
  R = numpy.empty(0)
  n = int(numpy.sqrt(Seeds)) + 4
  prev = []
  dev = x*dev
  while (len(S)<n**2):
    if (len(S)==0):
      S = numpy.append(S, numpy.array([[-2*x, -2*x]]), axis = 0)
      R = numpy.append(R, numpy.array([x/4]), axis = 0)
    a = (S[len(S)-1, :] + numpy.array([x, 0] + numpy.random.normal(0, dev, (1, 2))))
    b = (numpy.array([-2*x, S[len(S)-1, 1]]) + numpy.array([0, x] + numpy.random.normal(0, dev, (1, 2))))
    if S[len(S)-1, 0]< Side+2*x:
      if not prev:
        S0 = (S[len(S)-1, :] + numpy.array([x, 0] + numpy.random.normal(0, dev, (1, 2))))
      else:
        S0 = numpy.array([S[len(S)-1, 0], S[len(S)-prev[0], 1]])  + numpy.array([x, x] + numpy.random.normal(0, dev, (1, 2)))
    else:
      if not prev:
        prev.append(len(S)-1)
        S0 = (numpy.array([-2*x, -2*x]) + numpy.array([0, x] + numpy.random.normal(0, dev, (1, 2))))
      else:
        prev[0] = len(S)-1-prev[0]
        S0 = (numpy.array([-2*x, S[len(S)-prev[0], 1]]) + numpy.array([0, x] + numpy.random.normal(0, dev, (1, 2))))
      
    if numpy.all(numpy.sqrt(numpy.sum((S0-S)**2, axis = 1))-R- x/4>0):
      S = numpy.append(S, S0, axis = 0)
      R = numpy.append(R, numpy.array([x/4]), axis = 0)
  return S, R

def voronoi(Seeds, Side, dev, w):
  r = Side*numpy.sqrt(2/Seeds/numpy.sqrt(3))
  n2 = int(Side/r)+1
  n1 = int(Side*2/r/numpy.sqrt(3))+1
  R = numpy.ones((n1+8)*(n2+8))*0.5
  S = numpy.zeros(((n1+8)*(n2+8), 2))
  # x = 1
  # r = x*(0.5/numpy.tan(pi/12)+1/numpy.sqrt(3))
  w0 = r/4*numpy.array([1, 1, w, w])
  w1 = r/4*numpy.array([1, 1, w, w])
  w2 = r/4*numpy.array([1, 1, w, w])
  dev = r*dev
  a0 = [0, 0, -0.35, 0.35, -0.35, 0.35]
  a1 = [0, 0.5, 0.35, -0.35, 0.35, -0.35]
  a2 = [0, 0, -0.35,  0.35, -0.35, 0.35]
  b0 = [0, 0, 0.35, -0.35, 0.35, -0.35]
  b1 = [0, 0.5, -0.35, 0.35]
  b2 = [0, 0, 0.35, -0.35]
  deva = [1, 1, 0.2, 0.2]
  inda = numpy.zeros(((n2+8), (n1+8), 2))
  for i in (numpy.arange(0, n2+8, 3)):
      for j in numpy.arange(0, n1+8, 3):
        ind = 3
        if i+3<n2+8 and i+2<n2+8 and j+3<n1+8 and j+2<n2+8:
          inda[i, j, 0] = a0[ind]; inda[i, j+1, 0] = a1[ind]; inda[i, j+2, 0] =  a2[ind]
          inda[i+1, j, 0] = b0[ind]; inda[i+1, j+1, 0] = b1[ind]; inda[i+1, j+2, 0] =  b2[ind]
          inda[i+2, j, 0] = a0[ind]; inda[i+2, j+1, 0] = a1[ind]; inda[i+2, j+2, 0] =  a2[ind]
          inda[i:i+3, j:j+3, 1] = deva[ind]


  for i in range(n2+8):
      for j in range(n1+8):
        if ((j+1)%3 == 1) and (i+1)%3==1:
          ind = numpy.random.choice([0, 1, 2, 2])
        if (j+1)%3 == 1:
              S[(n1+8)*i+j, 0] = r*numpy.sin(pi/3)*j - 4*r*numpy.sin(pi/3)+numpy.random.uniform(-dev, dev)*inda[i, j, 1]
              a00 = b0[ind] if (i+1)%2==1 else a0[ind]
              S[(n1+8)*i+j, 1] = r*i - 4*r+numpy.random.uniform(-dev, dev) + inda[i, j, 0]*r
              R[(n1+8)*i+j] = w0[1]

        if (j+1)%3 == 2:
          S[(n1+8)*i+j, 0] = r*numpy.sin(pi/3)*j -4*r*numpy.sin(pi/3)+numpy.random.uniform(-dev, dev)*inda[i, j, 1]
          a10 = b1[ind] if (i+1)%2==1 else a1[ind]
          S[(n1+8)*i+j, 1] = inda[i, j, 0]*r + r*i - 4*r
          R[(n1+8)*i+j] = w1[1]
        if (j+1)%3 == 0:
              S[(n1+8)*i+j, 0] = r*numpy.sin(pi/3)*j - 4*r*numpy.sin(pi/3)+numpy.random.uniform(-dev, dev)*inda[i, j, 1]
              a20 = b2[ind] if (i+1)%2==1 else a2[ind]
              S[(n1+8)*i+j, 1] = r*i - 4*r+numpy.random.uniform(-dev, dev) + inda[i, j, 0]*r
              R[(n1+8)*i+j] = w2[1]
  return S, R
      # numpy.random.uniform(-dev, dev)