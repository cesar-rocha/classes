program pb3

    implicit none
    integer, parameter :: N = 512
    integer :: i, j
    real (kind=8), parameter :: pi = acos(-1.d0)
    real (kind=8) :: L, dx, dk, nu, dt
    real (kind=8), dimension(N) :: uini, u, x
    real (kind=8), dimension(N/2+1):: k, k2
    complex (kind=8), dimension(N/2+1) :: uhat, Luhat, NLuhat
    external :: NLBurgers

    ! open files to write to disk
    open(unit=20, file="Burgers.out", action="write", status="replace")


    ! Physical domain
    L = 2.d0*pi
    dx = L/N

    do i=0,N-1
        x(i+1) = i*dx
    enddo

    ! Fourier domain
    dk = 2.d0*pi/L          ! spectral resolution in rad / [L]
    call wavenumber(N,L,k)  ! wavenumber array in rad / [L]
    k2 = k**2

    ! Constant parameters
    nu = 1.d-3              ! viscosity coefficient [L]^2/[T]

    dt = 1.d-2

    ! initial condition
    do i=1,N
        uini(i) = sin(x(i))
    enddo

    call rfft(uini,uhat,N)

    ! march forward in time
    do j=1,100000 
        call stepforward(uhat,N/2+1,nu,k2,dt)
    enddo

    call irfft(uhat,u,N)


    ! write to disk
    do i=1,N
        write(20,*) x(i), uini(i), u(i)
        !write(*,*) real(uhat(i)),aimag(uhat(i))
    enddo


end program pb3


!------------------------------------------------------------
! subroutines
!------------------------------------------------------------


! function to evaluate the linear part
subroutine LBurgers(u,Lu,k2,nu,n)

    integer, intent(in) :: n
    real (kind=8), intent(in) :: nu
    real (kind=8), dimension(n), intent(in) :: k2
    complex (kind=8), dimension(n), intent(in) :: u
    complex (kind=8), dimension(n), intent(out) :: Lu

    Lu = -nu*k2*u

end subroutine LBurgers
  
! function to evaluate the linear part
subroutine NLBurgers(u,NLu,n)

    integer, intent(in) :: n
    complex (kind=8), dimension(n), intent(in) :: u
    complex (kind=8), dimension(n), intent(out) :: NLu

    NLu = 0.d0

end subroutine NLBurgers
 

! step forward
subroutine stepforward(u,n,nu,k2,dt)

    implicit none

    ! subroutine arguments
    integer, intent(in) :: n
    real (kind=8), intent(in) :: dt, nu
    real (kind=8), dimension(n), intent(in) :: k2
    complex (kind=8), dimension(n), intent(inout) :: u


    ! local variables
    real (kind=8), parameter :: a1 = 29.d0/96.d0, a2 = -3.d0/40.d0, &
                                    a3 = 1.d0/6.d0
    real (kind=8), parameter :: b1 = 37.d0/160.d0, b2 = 5.d0/24.d0, &
                                    b3 = 1.d0/6.d0
    real (kind=8), parameter :: c1 = 8.d0/15.d0, c2 = 5.d0/12.d0, &
                                    c3 = 3.d0/4.d0
    real (kind=8), parameter :: d1 = -17.d0/60.d0, d2 = 5.d0/12.d0
    real (kind=8), dimension(n) :: Lin
    complex (kind=8), dimension(n) :: NL1, NL2
    integer :: i


    interface
        subroutine  NLBurgers(u,NLu,n)
            integer, intent(in) :: n
            complex (kind=8), dimension(n), intent(in) :: u
            complex (kind=8), dimension(n), intent(out) :: NLu
        end subroutine NLBurgers
    end interface

    Lin = -nu*dt*k2

    ! step forward, computing updates in place
    call NLBurgers(u,NL1,n)
    u = ( (1.d0 + a1*Lin)/(1.d0 - b1*Lin) )*u + c1*dt*NL1

    !NL2 = NL1
    call NLBurgers(u,NL1,n)
    u = ( (1.d0 + a2*Lin)/(1.d0 - b2*Lin) )*u + c2*dt*NL1 + d1*dt*NL1

    !NL2 = NL1
    print *, maxval(real(NL2),n/2+1),  maxval(aimag(NL2),n/2+1)
    call NLBurgers(u,NL1,n)
    u = ( (1.d0 + a3*Lin)/(1.d0 - b3*Lin) )*u + c3*dt*NL1 + d2*dt*NL1

!    do i=1,n
!        u(i) = ( (1.d0 + a1*Lin(i))/(1.d0 - b1*Lin(i)) )*u(i) 
!        u(i) = ( (1.d0 + a2*Lin(i))/(1.d0 - b2*Lin(i)) )*u(i)
!        u(i) = ( (1.d0 + a3*Lin(i))/(1.d0 - b3*Lin(i)) )*u(i)
!    enddo
!

end subroutine stepforward

! wavenumber
subroutine wavenumber(N,L,k)

    implicit none
    integer, intent(in) :: N
    real (kind=8), intent(in) :: L
    real (kind=8), dimension(N/2+1), intent(out) :: k

    integer :: i

    do i=1,N/2+1
        k(i) = real(i-1)
    enddo

end subroutine wavenumber

! fft
subroutine rfft(x,xhat,n)

    ! It computes a 1D real to complex the forward 1D fft 
    !   out-of-place 
    !
    ! Input
    ! =======
    ! - x : is the array to be  transformed
    ! - n : the size of the array x
    !
    ! Output
    ! ========
    ! - xhat : the complex-valued DFT of x

    implicit none
    include 'fftw3.f'
    
    integer, intent(in) :: n
    real (kind=8), dimension(n), intent(in) :: x
    complex (kind=8), dimension(n/2+1), intent(inout) :: xhat

    integer (kind=8) :: plan

    ! compute forward 1d fft for real input
    call dfftw_plan_dft_r2c_1d(plan,n,x,xhat,FFTW_ESTIMATE)
    call dfftw_execute(plan)
    call dfftw_destroy_plan(plan)

end subroutine rfft

! ifft
subroutine irfft(xhat,x,n)

    ! It computes a 1D complex to real the forward 1D fft 
    !   out-of-place 
    !
    ! Input
    ! =======
    ! - xhat : is the array to be  transformed
    ! - n :    the size of the array x
    !
    ! Output
    ! ========
    ! - x : the iDFT of xhat

    implicit none
    include 'fftw3.f'
    
    integer, intent(in) :: n
    real (kind=8), dimension(n), intent(inout) :: x
    complex (kind=8), dimension(n/2+1), intent(in) :: xhat

    integer (kind=8) :: plan

    ! compute forward 1d fft for real input
    call dfftw_plan_dft_c2r_1d(plan,n,xhat,x,FFTW_ESTIMATE)
    call dfftw_execute(plan)
    call dfftw_destroy_plan(plan)

    x = x/n

end subroutine irfft





