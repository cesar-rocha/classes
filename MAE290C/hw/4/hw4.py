import numpy as np
from numpy import sqrt, sin ,cos ,pi
import matplotlib.pyplot as plt

def Thomas(a,b,c,s):
    ''' it solves T x = s where T is a tridiagonal matrix
            with constants coefficients

            Inputs
            ------------------------------------------
            a :: main diagonal coefficients
            b :: lower diagonals coefficients
            c :: upper diagonals coefficients
            s :: righ-hand side (input array)

            Outputs
            ------------------------------------------
            x :: solution array
    '''

    assert s.size == a.size, "s and a must have the same length"

    n = s.size

    rhs = s.copy()
    f = a.copy()                # assemble arrays: main diagonal
    e = b.copy()                #                  lower diagonal
    g = c.copy()                #                  upper diagonal
    x = np.empty(n)             #                  solution

    for i in range(1,n):
        factor = e[i-1]/f[i-1]   # perform Gaussian elimination
        f[i] -=  factor*g[i-1]
        rhs[i] -=  factor*rhs[i-1]
        del factor

    x[-1] = rhs[-1]/f[-1]   # backward substitution to solve for x
    for i in range(n-2,-1,-1):
        x[i] = (rhs[i]-g[i]*x[i+1])/f[i]

    return x

def Frhs(u,t,dt,gamma,A=.1,bc='a'):

    u[0] = sin(A*t)

    if bc == 'a':
        u[-1] = 0.
    elif bc == 'b':
        u[-1] = 2.*u[-2]-u[-3]
    elif bc == 'c':
        u[-1] = 3.*u[-2]-3.*u[-3]+u[-4]
    elif bc == 'd':
        u[-1] = u[-2]

    if bc == 'e':
        rhs = np.zeros(u.size-1)
        rhs[:-1] =  gamma*(u[:-2] - u[2:])
        rhs[-1] = 2.*gamma*(u[-2]-u[-1])
        rhs = np.hstack([u[0],rhs])
    elif bc == 'f':
        rhs = np.zeros(u.size-1)
        rhs[:-1] =  gamma*(u[:-2] - u[2:])
        rhs[-1] = gamma*(-u[-3] + 4.*u[-2]-3.*u[-1])
        rhs = np.hstack([u[0],rhs])
    else:
        rhs = gamma*(u[:-2] - u[2:])
        rhs = np.hstack([u[0],rhs,u[-1]])

    return rhs

def step_forward(u,t,dt,gamma,Frhs,A=.1,bc='a',\
                a1 = 8/15.,a2 = 2/3.):

    # the right-hand side
    #s = Frhs(u,t,dt,gamma,A=A,bc=bc)

    # the tridiagonal matrix
    #a = np.ones(u.size)
    #b = -gamma*np.ones(u.size-1)
    #c = -b.copy()
    #unew = Thomas(a,b,c,s)

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

    # physical constants
    c = 1.

    # time stepping
    dt = 0.01
    gamma = 0.5*c/dx
    tmax = 200.
    nplot = int(25/dt)
    ndt = 0
    t = 0.

    # initial condition
    u = np.zeros_like(x)
    A = 3.

    # boundary condition
    bc = 'f'

    # step forward
    while t<tmax:

        if (ndt%nplot==0):

            plt.figure(figsize=(8,6))
            plt.plot(x,u)
            plt.xlabel(r'$x$')
            plt.ylabel(r'$u$')
            plt.xlim(0,10.)
            plt.ylim(-1.,1.)
            plt.title('A = '+str(A)+', t = '+str(int(np.round(t))))
            figtit = bc + '_' + str(int(np.round(t)))
            plt.savefig('figs/3/'+figtit)
            plt.close(1)

        u = step_forward(u,t,dt,gamma,Frhs,A=A,bc=bc)
        t += dt
        ndt += 1
