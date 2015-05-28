program vort_eqn

    ! MAE290C, Winter 2015
    ! Final Project
    ! Cesar B. Rocha    

    use fft_mod

    implicit none
    integer :: nx, ny, i, j, nmax, nsave
    character (len=4):: fno
    real (kind=8), parameter :: pi = acos(-1.d0)
    real (kind=8) :: L, dx, dk, nu, nu4, dt, tmax, tsave, ke, t, cfl, dy
    real (kind=8), dimension(:,:), allocatable :: q, qi, p, kx, ky, kappa2, filt, x, y, kappa_nd
    complex (kind=8), dimension(:,:), allocatable :: qh, ph,qho, jh
    integer (kind=8):: plan_forward, plan_backward    
    logical :: use_filter

    ! flags
    use_filter = .true.

    ! parameters
    nx = 1024
    ny = 1024
    !nu = 1.d0/5.d4
    nu = 1.d-10
    !nu4 = 3.125d-8
    nu4 = 0.d0
    dt = 1.d-3
    dx = 2.d0*pi/nx
    dy = 2.d0*pi/ny

    print *, "dx = ", dx

    nmax = 50000
    nsave = 500

    ! allocate variable
    allocate(q(ny,nx), p(ny,nx),x(ny,nx),y(ny,nx))
    allocate(qh(ny/2+1,nx),ph(ny/2+1,nx),jh(ny/2+1,nx))    
    allocate(kx(ny/2+1,nx),ky(ny/2+1,nx),kappa2(ny/2+1,nx), filt(ny/2+1,nx), kappa_nd(ny/2+1,nx) )

    ! initialize fft stuff
    call init_plan_rfft2(ny,nx,plan_forward)
    call init_plan_irfft2(ny,nx,plan_backward)
    call rfft2_wavenumber(ny,nx,kx,ky)
    kappa2 = kx**2 + ky**2

    ! initial field
    ! a plane wave
    do i=1,nx
        x(:,i) = dx*(i-1.d0)
    end do
    do i=1,ny
        y(i,:) = dy*(i-1.d0)
    end do
    print *,maxval(x),maxval(y)

    call random_number(q) 
    !q = cos(4.d0*x + 4.d0*y)
    q = q-5.d-1
    q =  q - sum(q)/nx/ny
    call rfft2(q,qh,ny,nx,plan_forward)
    !qh(:,nx/3:2*nx/3) = 0.d0
    !qh(ny/3:,:) = 0.d0

    call invert(qh,ph,ny,nx,kappa2)    
    ph(1,1) = 0.d0
    call calc_diag(nx,ny,dx,dt,ph,kappa2,ke,cfl)

    print *, "KE = ", ke

    ph = ph/dsqrt(ke)

    call calc_diag(nx,ny,dx,dt,ph,kappa2,ke,cfl)    
    print *, "KE = ", ke
    print *, "Max allow dt = ", cfl*dx

    qh = -kappa2*ph     
    call invert(qh,ph,ny,nx,kappa2)    
    call calc_diag(nx,ny,dx,dt,ph,kappa2,ke,cfl)    
    print *, "KE = ", ke

    qi = q
 
    !print *, sum(q)/nx/ny

    ! set up filter or dealiasing
    
    if (use_filter) then
     kappa_nd = sqrt((kx*dx)**2.d0+(ky*dy)**2.d0)
     filt = exp(-23.6d0*(kappa_nd-0.65*pi)**4)  
     do i=1,nx
        do j=1,ny/2+1
            if (kappa_nd(j,i)<0.65*pi) then
                filt(j,i) = 1.d0
            end if
        end do
     end do

    else
        filt = 1.d0
        filt(:,nx/3:2*nx/3) = 0.d0  ! for dealiasing using 2/3 rule
        filt(ny/3:,:) = 0.d0
    end if

    open(unit=20, file="q", action="write", status="replace")
    open(unit=21, file="qi", action="write", status="replace")
    
    do j=1,nx
        write (20,*) q(j,:)
        write (21,*) qi(j,:)        
    enddo

   ! step forward in time
    do i=1,nmax

        call stepforward(ny,nx,kx,ky,kappa2,ph,qh,nu,nu4,dt,plan_forward,plan_backward,filt)
       
        ! print out diagnostics and save q to disk
        if (mod(i,nsave)==0.d0) then
            
            call calc_diag(nx,ny,dx,dt,ph,kappa2,ke,cfl)
            print *, "t = ", i*dt, "KE = ", ke, " CFL = ", cfl

            qho = qh
            call irfft2(qho,q,ny,nx,plan_backward)
            
            do j=1,nx
                write (20,*) q(j,:)
            enddo

        endif

    end do

end program vort_eqn

! from qh compute ph
subroutine invert(qh,ph,m,n,kappa2)

    implicit none
    integer, intent(in) :: m, n    
    complex (kind=8), dimension(m/2+1,n), intent(in) :: qh
    real (kind=8), dimension(m/2+1,n), intent(in) :: kappa2    
    complex (kind=8), dimension(m/2+1,n), intent(out) :: ph

    integer :: i

    ph = -qh/kappa2
    ph(1,1) = 0.d0

end subroutine invert

! calculate kinetic energy
subroutine calc_diag(m,n,dx,dt,ph,kappa2,ke,cfl)

    implicit none

    integer, intent(in) :: m, n 
    real (kind=8), intent(in) :: dt, dx
    complex (kind=8), dimension(m/2+1,n), intent(in) :: ph
    real (kind=8), dimension(m/2+1,n), intent(in) :: kappa2
    real (kind=8), intent(out) :: ke, cfl
   
    real (kind=8), dimension(m/2+1,n) :: ked 

    ked = real(kappa2*ph*conjg(ph))/m/m/n/n
    ke = sum(ked)
    cfl = maxval(dsqrt(ked))*dt/dx 
    
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
    complex (kind=8), dimension(m/2+1,n), intent(inout) :: ph, qh
    complex (kind=8), dimension(m/2+1,n), intent(out) :: jh
    real (kind=8), dimension(m/2+1,n), intent(in) :: k, l
    integer (kind=8), intent(in) :: planf, planb

    complex (kind=8), parameter :: ii = cmplx(0.d0, 1.d0) 
    complex (kind=8), dimension(m/2+1,n) :: uh, vh    
    real (kind=8), dimension(m,n) :: u, v, q

    ! compute velocities and vorticity in physical space
    
    ! for dealiasing
    !qh(:,n/3:2*n/3) = 0.d0
    !qh(m/3:,:) = 0.d0
    !ph(:,n/3:2*n/3) = 0.d0
    !ph(m/3:,:) = 0.d0

    uh = -ii*l*ph
    vh =  ii*k*ph

    call irfft2(uh,u,m,n,planb)
    call irfft2(vh,v,m,n,planb)
    call irfft2(qh,q,m,n,planb)

    call rfft2(u*q,uh,m,n,planf)
    call rfft2(v*q,vh,m,n,planf)

    ! now the product
    jh = -(ii*k*uh + ii*l*vh)

    !jh(:,n/3:2*n/3) = 0.d0
    !jh(m/3:,:) = 0.d0

end subroutine Jacobian 

! step forward
subroutine stepforward(m,n,k,l,kappa2,ph,qh,nu,nu4,dt,planf,planb,filt)

    implicit none

    ! subroutine arguments
    integer, intent(in) :: m, n
    real (kind=8), intent(in) :: dt, nu, nu4
    real (kind=8), dimension(m/2+1,n), intent(in) :: k, l, kappa2, filt
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
    real (kind=8), dimension(m/2+1,n) :: Lin
    complex (kind=8), dimension(m/2+1,n) :: jhold, jh, qho
    integer :: i

    interface
        subroutine  Jacobian(m,n,ph,qh,k,l,jh,planf,planb)
            integer, intent(in) :: m, n
            real (kind=8), dimension(m/2+1,n), intent(in) :: k, l
            complex (kind=8), dimension(m/2+1,n), intent(inout) :: ph, qh
            complex (kind=8), dimension(m/2+1,n), intent(out) :: jh
            integer (kind=8), intent(in) :: planf, planb            
        end subroutine Jacobian
    end interface

    interface
        subroutine  invert(qh,ph,m,n,kappa2)
            integer, intent(in) :: m, n    
            complex (kind=8), dimension(m/2+1,n), intent(in) :: qh
            real (kind=8), dimension(m/2+1,n), intent(in) :: kappa2    
            complex (kind=8), dimension(m/2+1,n), intent(out) :: ph
        end subroutine invert
    end interface

    Lin = -nu*dt*kappa2 - nu4*dt*(kappa2**2)
    !print *, Lin

    jh = 0.d0
    jhold = jh

    ! step forward, computing updates in place
    qho = qh
    call Jacobian(m,n,ph,qh,k,l,jh,planf,planb)
    qh = filt*( ( (1.d0 + a1*Lin)/(1.d0 - b1*Lin) )*qho + c1*dt*jh )

    jhold = jh
    qho = qh    
    call invert(qh,ph,m,n,kappa2) 
    call Jacobian(m,n,ph,qh,k,l,jh,planf,planb)
    qh = filt*( ( (1.d0 + a2*Lin)/(1.d0 - b2*Lin) )*qho + c2*dt*jh + d1*dt*jhold )

    jhold = jh
    qho = qh
    call invert(qh,ph,m,n,kappa2)
    call Jacobian(m,n,ph,qh,k,l,jh,planf,planb)
    qh = filt*( ( (1.d0 + a3*Lin)/(1.d0 - b3*Lin) )*qho + c3*dt*jh + d2*dt*jhold )

    call invert(qh,ph,m,n,kappa2)    

end subroutine stepforward

