program pb3

    ! MAE290C, Winter 2015
    ! Homework 2, Problem 3: 1D Burgers eqn.
    ! Cesar B. Rocha    

    implicit none
    integer :: N, i, j, nmax, nsave
    character (len=4):: fno
    real (kind=8), parameter :: pi = acos(-1.d0)
    real (kind=8) :: L, dx, dk, nu, nu6, dt, tmax,tsave
    real (kind=8), dimension(:), allocatable :: uini, u, x, k, filt
    complex (kind=8), dimension(:), allocatable :: uhat, uhats
    logical :: use_filter

    ! reading simulation params
    read (*,*) N            ! Number of points in physical space
    read (*,*) fno          ! string with N for output files
    read (*,*) nu           ! viscosity coefficient [L]^2/ [T] 
    read (*,*) nu6          ! 6th order hyperviscosity coefficient [L]^6/ [T] 
    read (*,*) use_filter   ! filter flag
    read (*,*) tmax         ! max time of simulation
    read (*,*) tsave        ! save every

    ! allocate variables based on N
    allocate(uini(N), u(N), x(N)) 
    allocate(k(N/2+1), uhat(N/2+1), uhats(N/2+1), filt(N/2+1))

    ! open files to write to disk

    open(unit=20, file="output/Burgers"//fno//".config", action="write", status="replace")
    open(unit=21, file="output/Burgers"//fno//".ini", action="write", status="replace")
    open(unit=22, file="output/Burgers"//fno//".phys", action="write", status="replace")
    open(unit=23, file="output/Burgers"//fno//".four", action="write", status="replace")
    open(unit=24, file="output/Burgers"//fno//".time", action="write", status="replace")
    open(unit=25, file="output/Burgers"//fno//".k", action="write", status="replace")

    ! Physical domain
    L = 2.d0*pi
    dx = L/N

    do i=0,N-1
        x(i+1) = i*dx
    enddo

    if (x(N)==L) then
        print *, "*** Warning: last element of x is redundant"
    endif

    ! Fourier domain
    dk = 2.d0*pi/L          ! spectral resolution in rad / [L]
    call wavenumber(N,L,k)  ! wavenumber array in rad / [L]

    ! spectral filter
    do i=1,N/2+1
        if (use_filter) then
            filt(i) = (1.d0+cos(2.d0*pi*(i-1)/N))/2.
        else
            filt(i) = 1.d0
        endif
    enddo

    ! Constant parameters
    !nu = 1.d-3              ! viscosity coefficient [L]^2/ [T]
    dt = .2*dx               ! dt_max: cfl < sqrt(3) for RK3
    nmax = ceiling(tmax/dt)
    nsave = ceiling(tsave/dt)
    write (20,*) "dx = ", dx
    write (20,*) "nu = ", nu
    write (20,*) "nu6 = ", nu6 
    write (20,*) "filter = ", use_filter 
    write (20,*) "tmax = ", tmax
    write (20,*) "dt = ", dt
    write (20,*) "nmax = ", nmax
    write (20,*) "nsave = ", nsave
    write (20,*) "kmax = ", k(N/2+1)

    ! initial condition
    do i=1,N
        uini(i) = sin(x(i))
        write(21,*) x(i),uini(i)   
    enddo

    write(25,*) k

    call rfft(uini,uhat,N)

    ! march forward in time
    do j=1,nmax

        ! filter
        uhat = filt*uhat

        if (mod(j,nsave)==0.d0) then
            ! some weird thing going on
            ! this doesn't work if I work directly on uhat
            ! i think fftw is changing uhat
            uhats = uhat
            call irfft(uhats,u,N)           
            print *, "t = ", j*dt, "CFL = ", maxval(abs(u))*dt/dx
            write (24,*) j*dt
            write (22,*) u
            write (23,*) dsqrt(real(uhat)**2 + aimag(uhat)**2 )
        endif

        call stepforward(uhat,N,nu,nu6,k,dt)

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
subroutine NLBurgers(uhat,k,NLuhat,n)

    ! Compute fully dealiased (3/2 rule) nonlinear term
    !  of 1D Burgers equation

    integer, intent(in) :: n
    complex (kind=8), dimension(n/2+1), intent(in) :: uhat
    complex (kind=8), dimension(n/2+1), intent(out) :: NLuhat
    real (kind=8), dimension(n/2+1), intent(in) :: k

    complex (kind=8), parameter :: ii = cmplx(0.d0, 1.d0) 
    complex (kind=8), dimension(3*n/4+1) :: uhat_pad, uxhat_pad, NLuhat_pad
    real (kind=8), dimension(3*n/2) :: u_pad, ux_pad

    interface
        subroutine  rfft(x,xhat,n)
            integer, intent(in) :: n
            real (kind=8), dimension(n), intent(in) :: x
            complex (kind=8), dimension(n/2+1), intent(inout) :: xhat
        end subroutine rfft
    end interface

    interface
        subroutine  irfft(xhat,x,n)
            integer, intent(in) :: n
            real (kind=8), dimension(n), intent(inout) :: x
            complex (kind=8), dimension(n/2+1), intent(in) :: xhat
        end subroutine irfft
    end interface

    ! initialize padded variables
    uhat_pad = 0.d0
    uxhat_pad = 0.d0
    ! pad u
    uhat_pad(1:n/2+1) = uhat
    call irfft(uhat_pad,u_pad,3*n/2)

    ! compute ux in Fourier
    uxhat_pad(1:n/2+1) = ii*k*uhat    
    call irfft(uxhat_pad,ux_pad,3*n/2)

    ! now compute the transform of the product
    call rfft(u_pad*ux_pad,NLuhat_pad,3*n/2)

     !NLuhat = 0.d0
    NLuhat = -NLuhat_pad(1:n/2+1)

end subroutine NLBurgers
 

! step forward
subroutine stepforward(u,n,nu,nu6,k,dt)

    implicit none

    ! subroutine arguments
    integer, intent(in) :: n
    real (kind=8), intent(in) :: dt, nu, nu6
    real (kind=8), dimension(n/2+1), intent(in) :: k
    complex (kind=8), dimension(n/2+1), intent(inout) :: u


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
        subroutine  NLBurgers(u,k,NLu,n)
            integer, intent(in) :: n
            real (kind=8), dimension(n), intent(in) :: k
            complex (kind=8), dimension(n), intent(in) :: u
            complex (kind=8), dimension(n), intent(out) :: NLu
        end subroutine NLBurgers
    end interface

    Lin = -nu*dt*(k**2) - nu6*dt*(k**6)

    ! step forward, computing updates in place
    call NLBurgers(u,k,NL1,n)
    u = ( (1.d0 + a1*Lin)/(1.d0 - b1*Lin) )*u + c1*dt*NL1

    NL2 = NL1
    call NLBurgers(u,k,NL1,n)
    u = ( (1.d0 + a2*Lin)/(1.d0 - b2*Lin) )*u + c2*dt*NL1 + d1*dt*NL2

    NL2 = NL1
    call NLBurgers(u,k,NL1,n)
    u = ( (1.d0 + a3*Lin)/(1.d0 - b3*Lin) )*u + c3*dt*NL1 + d2*dt*NL2

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





