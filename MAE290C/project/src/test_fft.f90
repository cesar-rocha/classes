program test

    ! This program showcases the use of the 2D real
    !  DFT from FFTW

    !use fft_mod, only: rfft2_ip, irfft2_ip, init_plan_rfft2
    use fft_mod

    implicit none

    integer :: n, i, j, nmax
    integer (kind=8):: plan_forward, plan_backward
    real (kind=8), allocatable :: A(:,:), Af(:,:)
    complex (kind=8), allocatable :: Ah(:,:)
    real (kind=8) :: k,t2,t1,t3,t4
    real (kind=8), parameter :: pi = acos(-1.d0)

    logical :: debug
    integer(kind=8) :: tclock1, tclock2, clock_rate

    n = 1024

    nmax = 1000

    allocate(A(n+2,n), Af(n+2,n))  ! need to pad two rows, since fortran stores
    allocate(Ah(n/2+1,n))          !    by column

    debug = .false.

    call random_number(A)

    ! keep old array to asses accuracy
    Af = A

    call cpu_time(t3)   ! start cpu timer

    ! initialize fft plans
    call init_plan_rfft2_ip(n,n,plan_forward)
    call init_plan_irfft2_ip(n,n,plan_backward)

    call cpu_time(t4)   ! end cpu timer


    call cpu_time(t1)   ! start cpu timer

    do j = 1,nmax

        !call rfft2(A,Ah,n,n,plan_forward)
        !call irfft2(Ah,A,n,n,plan_backward)

        call rfft2_ip(A,n,n,plan_forward)
        call irfft2_ip(A,n,n,plan_backward)

    enddo

    call cpu_time(t2)   ! end cpu timer


    print 9,  t4-t3
    9 format("Initialize plans in CPU time = ",f12.8, " seconds")

    print 10, nmax, t2-t1
    10 format("Performed ",i6.6," 2D FORWARD/BACKWARD RFFT2: CPU time = ",f12.8, " seconds")

    print 11, nmax, maxval(abs(A(1:n,:)-Af(1:n,:)))
    11 format("Absolute error after ",i6.6," 2D FORWARD/BACKWARD RFT2",1pe17.8)

    print 12, nmax, maxval(abs( (A(1:n,:)-Af(1:n,:))/Af(1:n,:) ))
    12 format("Relative error after ",i6.6," 2D FORWARD/BACKWARD RFT2",1pe17.8)


    open(unit=20, file="A", action="write", status="replace")
    open(unit=21, file="Af", action="write", status="replace")

    do j=1,n
        write (20,*) A(j,:)
        write (21,*) Af(j,:)
    enddo

end program test
