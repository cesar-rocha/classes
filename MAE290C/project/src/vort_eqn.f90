program vort_eqn

    ! MAE290C, Winter 2015
    ! Final Project
    ! Cesar B. Rocha    

    use fft_mod

    implicit none
    integer :: nx, ny, i, j, nmax, nsave, ntd
    character (len=4):: fno
    real (kind=8), parameter :: pi = acos(-1.d0)
    real (kind=8) :: L, dx, dk, nu, nu4, dt, tmax, tsave, ke, ens, t, cfl, dy
    real (kind=8), dimension(:,:), allocatable :: q, qi, p, kx, ky, kappa2,kappa2i
    real (kind=8), dimension(:,:), allocatable :: filt, x, y, kappa_nd, Lin
    real (kind=8), dimension(:,:), allocatable :: cmag, cphaser, cphasei
    complex (kind=8), dimension(:,:), allocatable :: qh, ph,qho, jh
    integer (kind=8):: plan_forward, plan_backward   
    complex (kind=8), parameter :: ii = cmplx(0.d0,1.d0)
    logical :: use_filter

    ! flags
    use_filter = .false.

    ! threads
    ntd = 4

    ! parameters
    nx =  1024
    ny = nx
    !nu = 0.d0*pi/1.d4
    nu = 0.d0
    nu = 2.d0*pi*2.d-5
    !nu4 = 3.125d-8
    nu4 = 0.d0
    dt = 1.d-3
    dx = 2.d0*pi/nx
    dy = 2.d0*pi/ny

    ! time params
    tmax = 2
    nmax = int(tmax/dt)
    nsave = 1000

    ! allocate variables
    allocate(q(ny,nx), p(ny,nx),x(ny,nx),y(ny,nx))
    allocate(qh(ny/2+1,nx),ph(ny/2+1,nx),jh(ny/2+1,nx))    
    allocate(kx(ny/2+1,nx),ky(ny/2+1,nx),kappa2(ny/2+1,nx),kappa2i(ny/2+1,nx))
    allocate(filt(ny/2+1,nx), kappa_nd(ny/2+1,nx), qi(ny,nx), Lin(ny/2+1,nx) )
    allocate(cmag(ny/2+1,nx),cphaser(ny/2+1,nx),cphasei(ny/2+1,nx))

    ! initialize fft stuff
    call init_plan_rfft2(ny,nx,ntd,plan_forward)
    call init_plan_irfft2(ny,nx,ntd,plan_backward)
    call rfft2_wavenumber(ny,nx,kx,ky)
    kappa2 = kx**2 + ky**2
    kappa2i = 1.d0/kappa2
    kappa2i(1,1) = 0.d0

    ! set up filter or dealiasing (should be a function)
    if (use_filter) then
     kappa_nd = dsqrt((kx*dx)**2.d0+(ky*dy)**2.d0)
     filt = dexp(-23.6d0*(kappa_nd-0.65d0*pi)**4)
    do i=1,nx
        do j=1,ny/2+1
            if (kappa_nd(j,i)<0.65d0*pi) then
                filt(j,i) = 1.d0
            end if
        end do
     end do
    else
        filt = 1.d0
        filt(:,nx/3:2*nx/3) = 0.d0  ! for dealiasing using 2/3 rule
        filt(ny/3:,:) = 0.d0
    end if
    filt(1,1) = 0.d0
    
    ! open units to write to disk
    open(unit=20, file="q", action="write", status="replace")
    open(unit=21, file="kappa2", action="write", status="replace")
    open(unit=22, file="qi", action="write", status="replace")
    open(unit=23, file="filt", action="write", status="replace")

    ! read initial condition 
    open(unit = 100, file = 'qi.txt', status = 'old', action = 'read')
     do j = 1,  nx
     read(100,*) qi(j,:)
    end do
    
    q = qi
    call rfft2(q,qh,ny,nx,plan_forward)

    do j=1,nx
        write (20,*) q(:,j)
        write (21,*) kappa2(:,j)
        write (22,*) qi(:,j)
        write (23,*) filt(:,j)
    enddo

    ! the linear operator 
    Lin = -1.d0*nu*dt*kappa2 -1.d0*dt*nu4*(kappa2**2)

    ! step forward in time
    do i=1,nmax

        call stepforward(ny,nx,kx,ky,kappa2,kappa2i,ph,qh,dt,plan_forward,&
                    plan_backward,filt,Lin)
       
        ! print out diagnostics and save q to disk
        if (mod(i,nsave)==0.d0) then
            
            call calc_diag(nx,ny,dx,dt,ph,kappa2,ke,ens)
            call calc_cfl(nx,ny,ph,dt,dx,kx,ky,cfl,plan_backward)
            
            print *, "t = ", i*dt, "KE = ", ke,"ENS = ", ens, " CFL = ", cfl

            qho = qh
            call irfft2(qho,q,ny,nx,plan_backward)
            
            do j=1,nx
                write (20,*) q(:,j)
            enddo

        endif

    end do

end program vort_eqn

! from qh compute ph
subroutine invert(qh,ph,m,n,kappa2i)

    implicit none
    integer, intent(in) :: m, n    
    complex (kind=8), dimension(m/2+1,n), intent(in) :: qh
    real (kind=8), dimension(m/2+1,n), intent(in) :: kappa2i   
    complex (kind=8), dimension(m/2+1,n), intent(out) :: ph

    integer :: i

    ph = -kappa2i*qh

end subroutine invert

! calculate kinetic energy
subroutine calc_diag(m,n,dx,dt,ph,kappa2,ke,ens)

    implicit none

    integer, intent(in) :: m, n 
    real (kind=8), intent(in) :: dt, dx
    complex (kind=8), dimension(m/2+1,n), intent(in) :: ph
    real (kind=8), dimension(m/2+1,n), intent(in) :: kappa2
    real (kind=8), intent(out) :: ke, ens
   
    real (kind=8), dimension(m/2+1,n) :: ked 

    ked = real(kappa2*ph*conjg(ph))/m/m/n/n
    ked(1,:) = ked(1,:)/2.d0
    ked(m/2+1,:) = ked(m/2+1,:)/2.d0
    ke = sum(ked)
    ens = sum(kappa2*ked)
    !cfl = maxval(dsqrt(ked))*dt/dx 
    
end subroutine calc_diag

! calculate CFL
subroutine calc_cfl(m,n,ph,dt,dx,k,l,cfl,planb)

    use fft_mod    

    implicit none

    integer, intent(in) :: m, n 
    real (kind=8), intent(in) :: dt, dx
    real (kind=8), intent(out) :: cfl
    complex (kind=8), dimension(m/2+1,n), intent(in) :: ph
    real (kind=8), dimension(m/2+1,n), intent(in) :: k, l
    integer (kind=8), intent(in) :: planb
    
    real (kind=8), dimension(m,n) :: u, v
    complex (kind=8), parameter :: ii = cmplx(0.d0, 1.d0) 
    
    call irfft2(-ii*l*ph,u,m,n,planb)
    call irfft2( ii*k*ph,v,m,n,planb)

    u = dsqrt(u**2 + v**2)
    cfl = maxval(u)*dt/dx

end subroutine calc_cfl

! function to evaluate the linear part
subroutine Jacobian(m,n,ph,qh,k,l,jh,planf,planb)

    ! Compute fully dealiased (2/3 rule) nonlinear term
    !  of 2D vort equation

    use fft_mod    

    implicit none

    integer, intent(in) :: m, n
    complex (kind=8), dimension(m/2+1,n), intent(in) :: ph, qh
    complex (kind=8), dimension(m/2+1,n), intent(out) :: jh
    real (kind=8), dimension(m/2+1,n), intent(in) :: k, l
    integer (kind=8), intent(in) :: planf, planb

    complex (kind=8), parameter :: ii = cmplx(0.d0, 1.d0) 
    complex (kind=8), dimension(m/2+1,n) :: uh, vh    
    real (kind=8), dimension(m,n) :: u, v, q

    ! compute velocities and vorticity in physical space
    uh = -ii*l*ph
    vh =  ii*k*ph

    call irfft2(uh,u,m,n,planb)
    call irfft2(vh,v,m,n,planb)
    call irfft2(qh,q,m,n,planb)

    call rfft2(u*q,uh,m,n,planf)
    call rfft2(v*q,vh,m,n,planf)

    ! the jacobian in fourier space
    jh = ii*k*uh + ii*l*vh

end subroutine Jacobian 

! step forward
subroutine stepforward(m,n,k,l,kappa2,kappa2i,ph,qh,dt,planf,planb,filt,Lin)

    implicit none

    ! subroutine arguments
    integer, intent(in) :: m, n
    real (kind=8), intent(in) :: dt
    real (kind=8), dimension(m/2+1,n), intent(in) :: k, l, kappa2, kappa2i, filt
    real (kind=8), dimension(m/2+1,n), intent(in) :: Lin
    complex (kind=8), dimension(m/2+1,n), intent(inout) :: ph, qh
    integer (kind=8), intent(in) :: planf, planb

    ! local variables
    real (kind=8), parameter :: a1 = 29.d0/96.d0, a2 = -3.d0/40.d0, &
                                    a3 = 1.d0/6.d0
    real (kind=8), parameter :: b1 = 37.d0/160.d0, b2 = 5.d0/24.d0, &
                                    b3 = 1.d0/6.d0
    real (kind=8), parameter :: c1 = 8.d0/15.d0, c2 = 5.d0/12.d0, &
                                    c3 = 3.d0/4.d0
    real (kind=8), parameter :: d1 = -17.d0/60.d0, d2 = -5.d0/12.d0
    complex (kind=8), dimension(m/2+1,n) :: jhold, jh, qho

    interface
        subroutine  Jacobian(m,n,ph,qh,k,l,jh,planf,planb)
            integer, intent(in) :: m, n
            real (kind=8), dimension(m/2+1,n), intent(in) :: k, l
            complex (kind=8), dimension(m/2+1,n), intent(in) :: ph, qh
            complex (kind=8), dimension(m/2+1,n), intent(out) :: jh
            integer (kind=8), intent(in) :: planf, planb            
        end subroutine Jacobian
    end interface

    interface
        subroutine  invert(qh,ph,m,n,kappa2i)
            integer, intent(in) :: m, n    
            complex (kind=8), dimension(m/2+1,n), intent(in) :: qh
            real (kind=8), dimension(m/2+1,n), intent(in) :: kappa2i  
            complex (kind=8), dimension(m/2+1,n), intent(out) :: ph
        end subroutine invert
    end interface

    ! step forward, computing updates in place
    qho = qh
    call invert(qho,ph,m,n,kappa2i)
    call Jacobian(m,n,ph,qh,k,l,jh,planf,planb)
    qh = filt*( ( (1.d0 + a1*Lin)/(1.d0 - b1*Lin) )*qho - c1*dt*jh )

    jhold = jh
    qho = qh    
    call invert(qho,ph,m,n,kappa2i) 
    call Jacobian(m,n,ph,qh,k,l,jh,planf,planb)
    qh = filt*( ( (1.d0 + a2*Lin)/(1.d0 - b2*Lin) )*qho - c2*dt*jh - d1*dt*jhold )

    jhold = jh
    qho = qh
    call invert(qho,ph,m,n,kappa2i)
    call Jacobian(m,n,ph,qh,k,l,jh,planf,planb)
    qh = filt*( ( (1.d0 + a3*Lin)/(1.d0 - b3*Lin) )*qho - c3*dt*jh - d2*dt*jhold )

end subroutine stepforward
