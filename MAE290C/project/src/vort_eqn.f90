program vort_eqn

    ! MAE290C, Winter 2015
    ! Final Project
    ! Cesar B. Rocha    

    use fft_mod

    implicit none
    integer :: nx, ny, i, j, nmax, nsave
    character (len=4):: fno
    real (kind=8), parameter :: pi = acos(-1.d0)
    real (kind=8) :: L, dx, dk, nu, nu6, dt, tmax, tsave
    real (kind=8), dimension(:,:), allocatable :: q, qi, p, kx, ky, kappa2, filt
    complex (kind=8), dimension(:,:), allocatable :: qh, ph
    integer (kind=8):: plan_forward, plan_backward    
    logical :: use_filter

    ! flags
    use_filter = .false.

    ! parameters
    nx = 32
    ny = 32
    nu = 1.d-4
    nu6 = 0.d0
    dt = 0.d-5

    nmax = 10

    ! allocate variable
    allocate(q(ny,nx), p(ny,nx))
    allocate(qh(ny/2+1,nx),ph(ny/2+1,nx))    
    allocate(kx(ny/2+1,nx),ky(ny/2+1,nx),kappa2(ny/2+1,nx))

    ! initialize fft stuff
    call init_plan_rfft2_ip(ny,nx,plan_forward)
    call init_plan_irfft2_ip(ny,nx,plan_backward)
    call rfft2_wavenumber(ny,nx,kx,ky)  
    kappa2 = kx**2 + ky**2

    ! initial field
    call random_number(q)
    qi = q
    call rfft2(q,qh,ny,nx,plan_forward)
    call invert(qh,ph,ny,nx,kappa2)

    ! step forward in time
    do i=1,1
        call stepforward(ny,nx,kx,ky,kappa2,ph,qh,nu,nu6,dt,plan_forward,plan_backward)
    end do

    !print *,qh
    call irfft2(qh,q,ny,nx,plan_backward)
    !print *, q-qi


end program vort_eqn


! from qh compute ph
subroutine invert(qh,ph,m,n,kappa2)

    implicit none
    integer, intent(in) :: m, n    
    complex (kind=8), dimension(m/2+1,n), intent(in) :: qh
    real (kind=8), dimension(m/2+1,n), intent(in) :: kappa2    
    complex (kind=8), dimension(m/2+1,n), intent(out) :: ph

    integer :: i

    do i=1,m/2+1
        if (i==1) then
            ph(i,:) = -qh(i,2:)/kappa2(i,2:)
            ph(i,i) = 0.d0
        else
            ph(i,:) = -qh(i,:)/kappa2(i,:)
        endif
    end do

end subroutine invert

! function to evaluate the linear part
subroutine Jacobian(m,n,ph,qh,k,l,jh,planf,planb)

    ! Compute fully dealiased (2/3 rule) nonlinear term
    !  of 2D vort equation ** NOT DEALIASED YET

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
    call irfft2(-ii*l*ph,u,m,n,planb)
    call irfft2(ii*k*ph,v,m,n,planb)
    call irfft2(qh,q,m,n,planb)

    u = u*q
    v = v*q

    call rfft2(u,uh,m,n,planf)
    call rfft2(v,vh,m,n,planf)

    ! now the product
    jh = -(ii*k*uh + ii*l*vh)*0.d0

end subroutine Jacobian 

! step forward
subroutine stepforward(m,n,k,l,kappa2,ph,qh,nu,nu6,dt,planf,planb)

    implicit none

    ! subroutine arguments
    integer, intent(in) :: m, n
    real (kind=8), intent(in) :: dt, nu, nu6
    real (kind=8), dimension(m/2+1,n), intent(in) :: k, l, kappa2
    complex (kind=8), dimension(m/2+1,n), intent(inout) :: ph, qh
    integer (kind=8), intent(in) :: planf, planb


    ! local variables
    real (kind=8), parameter :: a1 = 29.d0/96.d0, a2 = -3.d0/40.d0, &
                                    a3 = 1.d0/6.d0
    real (kind=8), parameter :: b1 = 37.d0/160.d0, b2 = 5.d0/24.d0, &
                                    b3 = 1.d0/6.d0
    real (kind=8), parameter :: c1 = 8.d0/15.d0, c2 = 5.d0/12.d0, &
                                    c3 = 3.d0/4.d0
    real (kind=8), parameter :: d1 = -17.d0/60.d0, d2 = 5.d0/12.d0
    real (kind=8), dimension(m/2+1,n) :: Lin
    complex (kind=8), dimension(m/2+1,n) :: jhold, jh
    integer :: i

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
        subroutine  invert(qh,ph,m,n,kappa2)
            integer, intent(in) :: m, n    
            complex (kind=8), dimension(m/2+1,n), intent(in) :: qh
            real (kind=8), dimension(m/2+1,n), intent(in) :: kappa2    
            complex (kind=8), dimension(m/2+1,n), intent(out) :: ph           
        end subroutine invert
    end interface

    Lin = -nu*dt*(kappa2) - nu6*dt*(kappa2**3)
    print *, Lin

    ! step forward, computing updates in place
    call Jacobian(m,n,ph,qh,k,l,jh,planf,planb)
    qh = ( (1.d0 + a1*Lin)/(1.d0 - b1*Lin) )*qh + c1*dt*jh

    jhold = jh
    call invert(qh,ph,m,n,kappa2)    
    call Jacobian(m,n,ph,qh,k,l,jh,planf,planb)
    qh = ( (1.d0 + a2*Lin)/(1.d0 - b2*Lin) )*qh + c2*dt*jhold + d1*dt*jh

    jhold = jh
    call invert(qh,ph,m,n,kappa2)    
    call Jacobian(m,n,ph,qh,k,l,jh,planf,planb)
    qh = ( (1.d0 + a3*Lin)/(1.d0 - b3*Lin) )*qh + c3*dt*jhold + d2*dt*jh

    call invert(qh,ph,m,n,kappa2)    



end subroutine stepforward

