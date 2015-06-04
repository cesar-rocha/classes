# MAE 290C, HW 4
# Cesar B Rocha

# This code solves the 1d wave equation
#   with various artificial computational BC

import numpy as np
from numpy import sqrt, sin ,cos ,pi
import matplotlib.pyplot as plt

def Frhs(u,t,dt,gamma,A=.1,bc=1):

    u[0] = sin(A*t)

    if bc == 1:
        u[-1] = 0.
    elif bc == 2:
        u[-1] = 2.*u[-2]-u[-3]
    elif bc == 3:
        u[-1] = 3.*u[-2]-3.*u[-3]+u[-4]
    elif bc == 4:
        u[-1] = u[-2]

    if bc == 5:
        rhs = np.zeros(u.size-1)
        rhs[:-1] =  gamma*(u[:-2] - u[2:])
        rhs[-1] = 2.*gamma*(u[-2]-u[-1])
        rhs = np.hstack([u[0],rhs])
    elif bc == 6:
        rhs = np.zeros(u.size-1)
        rhs[:-1] =  gamma*(u[:-2] - u[2:])
        rhs[-1] = gamma*(-u[-3] + 4.*u[-2]-3.*u[-1])
        rhs = np.hstack([u[0],rhs])
    else:
        rhs = gamma*(u[:-2] - u[2:])
        rhs = np.hstack([u[0],rhs,u[-1]])

    return rhs

def step_forward(u,t,dt,gamma,Frhs,A=.1,bc=1,\
                a1 = 8/15.,a2 = 2/3.):

    f1 = Frhs(u,t,dt,gamma,A=A,bc=bc)
    u1 = u + dt*a1*f1
    f2 = Frhs(u1,t,dt,gamma,A=A,bc=bc)
    u2 = u + dt*a2*f2
    f3 = Frhs(u2,t,dt,gamma,A=A,bc=bc)
    unew = u + dt/4.*(f1+3.*f3)

    return unew

if __name__ == "__main__":

    plt.close('all')

    # physical domain
    N = 200
    dx = 10./N
    x = np.arange(0,10.+dx,dx) # including bcs

    # physical bc frequency
    A = 3.

    # physical constants
    c = 1.

    # time stepping
    dt = 0.01
    gamma = 0.5*c/dx
    tmax = 35.
    nplot = int(5/dt)
    ndt = 0
    t = 0.

    # boundary condition
    for bc in range(7):

        t, ntd = 0.,0
        u = np.zeros_like(x)
        
        # step forward
        while t<=tmax:

            if (ndt%nplot==0):

                plt.figure(figsize=(8,6))
                plt.plot(x,u)
                plt.xlabel(r'$x$')
                plt.ylabel(r'$u$')
                plt.xlim(0,10.)
                plt.ylim(-1.,1.)
                plt.title('BC (' + str(bc+1) + '),' + ' A = ' + str(A) + ', c = '+ str(c)+\
                        ', t = ' +str(int(np.round(t))))
                figtit = str(bc+1) + '_' + str(int(np.round(t)))
                plt.savefig('figs/3/'+figtit)
                plt.close(1)

            u = step_forward(u,t,dt,gamma,Frhs,A=A,bc=bc+1)
            t += dt
            ndt += 1
