program pb1

    implicit none
    real (kind=8) :: u0,dt,pi,tmax,t
    complex (kind=16) :: lam,lamdt,u
    integer :: nmax,n


    open(unit=12, file="aoutput.txt", action="write", status="replace")


    ! pi
    pi = acos(-1.d0)

    ! time
    dt = .5
    tmax = 10.d0
    nmax = ceiling(tmax/dt)

    print 11, nmax
    11 format('nmax=', i10)

    u0 = 1.d0
    u = cmplx(u0,0.d0)

    t = 0.d0

    lamdt = cmplx(pi/2.d0,0.5d0)

    write(12,*) t, real(u),aimag(u)

    do n=1,nmax
        call stepforward(u,lamdt)
        t = t+dt
        write(12,*) t,real(u),aimag(u)
    enddo


end program pb1

!------------------------------------------------------------
! subroutine to step forward using RKW3
!------------------------------------------------------------

subroutine stepforward(u,lamdt)

    implicit none

    ! subroutine arguments
    complex (kind=16), intent(in) :: lamdt
    complex (kind=16), intent(inout) :: u

    ! local variables
    real (kind=8), parameter :: a1 = 29.d0/96.d0, a2 = -3.d0/40.d0, & 
                                    a3 = 1.d0/6.d0
    real (kind=8), parameter :: b1 = 37.d0/160.d0, b2 = 5.d0/24.d0, & 
                                    b3 = 1.d0/6.d0

    ! step forward, computing updates in place
    u = ( (1.d0-lamdt*a1)/(1.d0 + lamdt*b1) )*u 
    u = ( (1.d0-lamdt*a2)/(1.d0 + lamdt*b2) )*u 
    u = ( (1.d0-lamdt*a3)/(1.d0 + lamdt*b3) )*u 

end subroutine stepforward

