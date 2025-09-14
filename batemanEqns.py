## Bateman Equations
# ------------------------------------------------------------------------------
# Adam Buchalter
#
# Analytical Solution to the Bateman Equations for a uniform decay constant and
# a stable terminal isotope
# ------------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
import math
import sys

# myLambda = 10*np.log(2)
myLambda = 1
N = 5

A = myLambda*np.eye(N, N, -1) - myLambda*np.eye(N, N, 0)
A[N-1,N-1]=0

# print(A)
# print(A@A@A@A@A@A)

# sys.exit("Error message")

# def N1(t, myLambda):
#     return np.exp(-myLambda*t)

# def N2(t, myLambda):
#     return myLambda*t*np.exp(-myLambda*t)

# def N3(t, myLambda):
#     return np.power(myLambda,2)*np.power(t,2)*np.exp(-myLambda*t)/2

def Nk(t, myLambda, k):
    return np.power(myLambda,k-1)*np.power(t,k-1)*np.exp(-myLambda*t)/math.factorial(k-1)

def Nfinal(t, myLambda, k):
    coeff = 0
    for n in range(k+1):
        coeff -= np.power(myLambda,n)*np.power(t,n)/math.factorial(n)
        print(n)
        print(np.power(myLambda,n)/math.factorial(n))
    return 1 + coeff*np.exp(-myLambda*t)

def myt(myLambda, k):
    return np.power(math.factorial(k-1),1/(k-1))/myLambda

N1_0 = 1
myTime = np.linspace(0, 10, num=1000)

for myk in range(1,N+1):
    plt.plot(myTime, N1_0*Nk(myTime, myLambda, myk))
plt.title(str(N)+' Unstable Isotopes')
plt.xlabel('Time [half-lives]')
plt.ylabel('Relative Population')
plt.savefig('results//bateman1.png')
plt.show()
plt.close()

for myk in range(1,N):
    plt.plot(myTime, N1_0*Nk(myTime, myLambda, myk))
plt.plot(myTime, N1_0*Nfinal(myTime, myLambda, N))
plt.title(str(N-1)+' Unstable Isotopes with a Stable Terminal Isotope')
plt.xlabel('Time [half-lives]')
plt.ylabel('Relative Population')
plt.savefig('results//bateman2.png')
plt.show()
plt.close()