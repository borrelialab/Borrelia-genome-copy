#!/usr/bin/env python
# mutatewlag simulates our infection by evolving the switch rates and including the four generation lag

import random

def mutate(n,m,k,lvariable,P):
    """
        Our function "mutate" will start with a given infection, P, and corresponding matrix of switch rates, lvariable.
        It will randomly mutate one of the matrix elements into a number in the range 9*10^-3 to 1*10^-6.
        It will simulate an infection with this matrix.  If this new infection is more fit, we keep the change.
        Otherwise, we change our matrix back. The function outputs the matrix we end up with and the infection, including the individual serotypes.  
    """
    summation1=0
    summation2=0
    for i in xrange(m):
        summation1=summation1+P[i]    
   
    A=np.ones(m)
    for i in xrange(m):
        P[i]=0
        A[i]=0
    varone=random.randint(0,n**2-1)
    mant=random.randint(mantBaseMin,mantBaseMax)
    exp=random.randint(expMin,expMax)
    save=lvariable[varone]
    lvariable[varone]=mant*10**(-exp)
    
    s=(n,m)
    p = np.ones(s)
    for i in xrange(n):
        p[i,0]=0
    p[0]=1

    a2=np.ones(n)
    for i in xrange(n):
        a2[i]=0	
    a = np.ones(s)
    
    for i in xrange(n):
        for j in range(m):
            a[i,j]=0
        
   
    for i in xrange(m-1):
       for j in xrange(n):
           x=0
           y=0
           for l in xrange(n):
              if j!=l:
                   if i <= mutateGeneration/k-1:
                       x=0
                       y=0
                   if i > mutateGeneration/k-1:
                       x=x+(p[j,i-mutateGeneration/k]*lvariable[j*n+l]) #the product of the amount of cells of serotype j that were present four generations prior, with the switch rate from j to l, summed over every l.  
                       y=y+(p[l,i-mutateGeneration/k]*lvariable[l*n+j]) #the product of the amount of cells of serotype l that were present four generations prior, with the switch rate from l to j, summed over every l.  
                       
           # We use a finite difference method with Forward Euler time stepping to solve our system of equations. The antibody response to a serotype will not be triggered until there are 1000 cells of that serotype.

          
           slope=k*(p[j,i]*(r*(1-P[i]/c)-kappa*abeffect[j,i])-x+y)
           if p[j,i]+slope>=0:
               p[j,i+1]=p[j,i]+slope
           if p[j,i]+slope<0:
               p[j,i+1]=0
        
           a[j,i+1]=a[j,i]+k*(a[j,i]*(rho*p[j,i]/(p[j,i]+phi)-u)+pi)
        
           abeffect[j,i+1]=a[j,i+1]*(1+a[j,i+1]/beta)
        
           P[i+1]=P[i+1]+p[j,i+1]    
    for i in xrange(m):
        summation2=P[i]+summation2
        
    if summation1>summation2:
        lvariable[varone]=save
    return [lvariable,p,P]
	

# we must initialize our program by simulating a reaction with random switch rates
import numpy as np
import pylab as pl

print "Please enter the number of serotypes [12]:"
n = raw_input(">") #12, number of serotypes
if n == "":
    n = 12
else:
    n = int(n)
    
print "Please enter the number of time steps [300]:"
m = raw_input(">") #300, number of time steps
if m == "":
    m = 300
else:
    m = int(m)
    
print "Please enter the length of a time step [0.1]:"
k = raw_input(">")
if k == "":
    k = 0.1
else:
    k = float(k)
    
print "Please enter the antibody synergy factor [0.5]:" 
beta = raw_input(">")
if beta == "":
    beta = 0.5
else:
    beta = float(beta)
    
print "Please enter the minimum base value for the mutation rate [1]:" 
mantBaseMin = raw_input(">")
if mantBaseMin == "":
    mantBaseMin = 1
else:
    mantBaseMin = int(mantBaseMin)
    
print "Please enter the maximum base value for the mutation rate [9]:" 
mantBaseMax = raw_input(">")
if mantBaseMax == "":
    mantBaseMax = 9
else:
    mantBaseMax = int(mantBaseMax)

print "Please enter the minimum negative exponent value for the mutation rate [3]:" 
expMin = raw_input(">")
if expMin == "":
    expMin = 3
else:
    expMin = int(expMin)
    if expMin < 0:
        expMin =  -1 * expMin
        print "Please remember that this is already assumed negative. Using ", expMin
    
print "Please enter the maximum negative exponent value for the mutation rate [6]:" 
expMax = raw_input(">")
if expMax == "":
    expMax = 6
else:
    expMax = int(expMax)
    if expMax < 0:
        expMax =  -1 * expMax
        print "Please remember that this is already assumed negative. Using ", expMax

print "Please enter the number for how many generations the infection takes to mutate [4]:" 
mutateGeneration = raw_input(">")
if mutateGeneration == "":
    mutateGeneration = 4
else:
    mutateGeneration = int(mutateGeneration)
    if mutateGeneration < 0:
        mutateGeneration =  -1 * mutateGeneration
        print "This number can not be negative. Using ", mutateGeneration

print "Please enter the ratio for the amount of antibody to serotype, at first generation [0.001]:" 
pi = raw_input(">")
if pi == "":
    pi = 0.001
else:
    pi = float(pi)
    if pi < 0:
        pi =  -1 * pi
        print "This number can not be negative. Using ", pi

print "Please enter the amount of cells of a certain serotype that would evoke half of the maximal immune response [52]:" 
phi = raw_input(">")
if phi == "":
    phi = 52 #5*10^6
else:
    phi = int(phi)
    if phi < 0:
        phi =  -1 * phi
        print "This number can not be negative. Using ", phi

print "Please enter the rate of loss of antibody due to degradation [0.001]:" 
u = raw_input(">")
if u == "":
    u = 0.001
else:
    u = float(u)
    if u < 0:
        u =  -1 * u
        print "This number can not be negative. Using ", u

print "Please enter the maximum cell density the host can support [10000000]:" 
c = raw_input(">")
if c == "":
    c = 10**7
else:
    c = int(c)
    if c < 0:
        c =  -1 * c
        print "This number can not be negative. Using ", c

print "Please enter the death rate constant [0.0001]:" 
kappa = raw_input(">")
if kappa == "":
    kappa = 0.0001
else:
    kappa = float(kappa)
    if kappa < 0:
        kappa =  -1 * kappa
        print "This number can not be negative. Using ", kappa

print "Please enter the undamped growth rate of bacteria [2]:" 
r = raw_input(">")
if r == "":
    r = 2
else:
    r = int(r)
    if r < 0:
        r =  -1 * r
        print "This number can not be negative. Using ", r

print "Please enter the relative growth rate of antibody response [1.6]:" 
rho = raw_input(">")
if rho == "":
    rho = 1.6 
else:
    rho = float(rho)
    if rho < 0:
        rho =  -1 * rho
        print "This number can not be negative. Using ", rho

print "Please enter the number of times you would like to run mutations [10000]:" 
mutationRuns = raw_input(">")
if mutationRuns == "":
    mutationRuns = 10000 
else:
    mutationRuns = int(mutationRuns)
    if mutationRuns < 0:
        mutationRuns =  -1 * mutationRuns
        print "This number can not be negative. Using ", mutationRuns
        
print "Please enter the length of time for one generation, in hours [8]:" 
genTime = raw_input(">")
if genTime == "":
    genTime = 8 
else:
    genTime = int(genTime)
    if genTime < 0:
        genTime =  -1 * genTime
        print "This number can not be negative. Using ", genTime
        
print "Running program with..."
print n,m,k,beta,mantBaseMin,mantBaseMax,expMin,expMax,mutateGeneration,pi,phi,u,c,kappa,r,rho,mutationRuns,genTime

P=np.ones(m) # total number of B. hermsii cells
A=np.ones(m) # total antibody reaction
for i in xrange(m):
    P[i]=0
    A[i]=0

lvariable= np.ones(n**2) #matrix of switch rates
for i in xrange(n**2):
    lvariable[i]=0

for i in xrange(n):
    for j in xrange(n):
        mant=random.randint(mantBaseMin,mantBaseMax)
	exp=random.randint(expMin,expMax)
	lvariable[j*n+i]=mant*10**(-exp)

s = (n,m)
p = np.ones(s) # matrix where each row corresponds to the infections of a certain serotype
for i in xrange(n):
    p[i,0]=0
p[0,0]=1 

a2=np.ones(n)
for i in xrange(n):
	a2[i]=0

a = np.ones(s) # matrix where each row corresponds to the antibody response to a certain serotype
abeffect=np.ones(s)

for i in xrange(n):
    for j in xrange(m):
        a[i,j]=0
        abeffect[i,j]=0

b2=np.ones(n)
for i in range(n):
	b2[i]=0

for i in xrange(m-1):    
    for j in xrange(n):
        x=0
        y=0        
        for l in xrange(n):
            if j!=l:
                if i <= mutateGeneration/k-1:
                    x=0
                    y=0
                if i > mutateGeneration/k-1:
                    x=x+(p[j,i-mutateGeneration/k]*lvariable[j*n+l])
                    y=y+(p[l,i-mutateGeneration/k]*lvariable[l*n+j])  # switch rates demonstrate lag #time of 4 generations        
        
        slope=k*(p[j,i]*(r*(1-P[i]/c)-kappa*abeffect[j,i])-x+y)
        
        if p[j,i]+slope>=0:
            p[j,i+1]=p[j,i]+slope
        
        if p[j,i]+slope<0:
            p[j,i+1]=0
                  
        a[j,i+1]=a[j,i]+k*(a[j,i]*(rho*p[j,i]/(p[j,i]+phi)-u)+pi)        
        abeffect[j,i+1]=a[j,i+1]*(1+a[j,i+1]/beta)        
        P[i+1]=P[i+1]+p[j,i+1]

# We run "mutate" many times to produce a fit matrix of switch rates.
for i in xrange(mutationRuns):
    [lvariable,p,P]=mutate(n,m,k,lvariable,P)
    
x=np.ones(m)
for i in range(m):
    x[i]=i*k*genTime #Set each generation's duration, in hours.
    
z=np.ones(m)
for i in range(m):
    z[i]=1000

# Calculate the log for each Y coordinate
P3=np.ones(m)
for i in range(m):
    P3[i]=np.log(P[i])

# Create CSV file and write headers
import csv
csvFile = csv.writer(open("MutateWLag_Plot_Data.csv", "wb"))
csvFile.writerow(["X","Y"])

# Clear new plot, plot points, and save X,Y coordinates to a CSV file
pl.clf()
for i in range(n):    
    pl.plot(x,p[i,0:])
    csvFile.writerow([x,p[i,0:]])    
    #pl.plot(x,z)    
    pl.plot(x,P)
    csvFile.writerow([x,P])    
    #pl.plot(x,a[i,0:])    
    #pl.plot(x,A)    
# Save plotted points as graph, in a postscript file
pl.savefig('MutateWLag_Graph.ps')

# Clear new plot, plot points, and save X,Y coordinates to a CSV file
pl.clf()
pl.semilogy(x, P3)
#csvFile.writerow([x, P3[i]])
# Save plotted points as graph, in a postscript file
pl.savefig('MutateWLag_Graph_Ylog.ps')