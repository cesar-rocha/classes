! $MAE290C/project/src/fft_mod.f90

module fft_mod

    use omp_lib    

    implicit none

    include 'fftw3.f'

    contains

!
!  Plan initializations
!

  ! plan fft2
  subroutine init_plan_rfft2_ip(m,n,ntd,plan)

    implicit none

    integer, intent(in) :: m, n, ntd
    integer (kind=8), intent(inout) :: plan

    real (kind=8), allocatable :: x(:,:)
    integer :: iret

    allocate(x(n+2,n))
  
    call dfftw_init_threads(iret)
    call dfftw_plan_with_nthreads(ntd)

    ! plan the forward 2d fft for real input, in place
    call dfftw_plan_dft_r2c_2d(plan,m,n,x,x,FFTW_PATIENT)

    deallocate(x)

  end subroutine init_plan_rfft2_ip

  ! plan ifft2
  subroutine init_plan_irfft2_ip(m,n,ntd,plan)

    implicit none

    integer, intent(in) :: m, n, ntd
    integer (kind=8), intent(inout) :: plan

    real (kind=8), allocatable :: xhat(:,:)
    integer :: iret
      
    allocate(xhat(n+2,n))

    call dfftw_init_threads(iret)
    call dfftw_plan_with_nthreads(ntd)

    ! plan the backward 2d fft for real input, in place
    call dfftw_plan_dft_c2r_2d(plan,m,n,xhat,xhat,FFTW_PATIENT)

    deallocate(xhat)

  end subroutine init_plan_irfft2_ip

! plan fft2
subroutine init_plan_rfft2(m,n,ntd,plan)

    implicit none

    integer, intent(in) :: m, n, ntd
    integer (kind=8), intent(inout) :: plan

    real (kind=8), allocatable :: x(:,:)
    complex (kind=8), allocatable :: xhat(:,:)

    integer :: iret
      
    allocate(x(m,n),xhat(m/2+1,n))

    call dfftw_init_threads(iret)
    call dfftw_plan_with_nthreads(ntd)

    ! plan forward 2d fft for real input
    call dfftw_plan_dft_r2c_2d(plan,m,n,x,xhat,FFTW_PATIENT)
    
    deallocate(x,xhat)

end subroutine init_plan_rfft2

! plan ifft2_ip
subroutine init_plan_irfft2(m,n,ntd,plan)

    implicit none

    integer, intent(in) :: m, n, ntd
    integer (kind=8), intent(inout) :: plan

    real (kind=8), allocatable :: x(:,:)
    complex (kind=8), allocatable :: xhat(:,:)
    
    integer :: iret
      
    allocate(x(m,n),xhat(m/2+1,n))

    call dfftw_init_threads(iret)
    call dfftw_plan_with_nthreads(ntd)

    ! plan forward 2d fft for real input
    call dfftw_plan_dft_c2r_2d(plan,m,n,xhat,x,FFTW_PATIENT)

    deallocate(x,xhat)

end subroutine init_plan_irfft2

!
!   2D DFTs
!

    ! fft2
    subroutine rfft2(x,xhat,m,n,plan)

        ! It computes a 2D real to complex the forward 2D fft
        !   out-of-place
        !
        ! Input
        ! =======
        ! - x : is the array to be  transformed
        ! - n : the size of the array x
        ! - plan: plan of execution initialized by init_plan_fft2
        ! Output
        ! ========
        ! - xhat : the complex-valued DFT of x

        implicit none

        integer, intent(in) :: m, n
        real (kind=8), dimension(m,n), intent(in) :: x
        complex (kind=8), dimension(m/2+1,n), intent(out) :: xhat
        integer (kind=8), intent(in) :: plan

        ! compute forward 2d fft for real input
        call dfftw_execute_dft_r2c(plan,x,xhat)

    end subroutine rfft2

    ! ifft2
    subroutine irfft2(xhat,x,m,n,plan)

        ! It computes a 2D complex to real the forward 2D fft
        !   out-of-place
        !
        ! Input
        ! =======
        ! - xhat : is the array to be  transformed
        ! - n :    the size of the array x
        ! - plan: plan of execution initialized by init_plan_irfft
        ! Output
        ! ========
        ! - x : the iDFT of xhat

        implicit none

        integer, intent(in) :: m,n
        real (kind=8), dimension(m,n), intent(out) :: x
        complex (kind=8), dimension(m/2+1,n), intent(in) :: xhat
        integer (kind=8), intent(in) :: plan

        ! compute forward 2d fft for real input
        call dfftw_execute_dft_c2r(plan,xhat,x)
        x = x/n/n

    end subroutine irfft2

    ! fft2 in place
    subroutine rfft2_ip(x,m,n,plan)

        ! It computes a 2D real to complex the forward 2D fft
        !   in-place. The output overwrite the input, and is
        !   therefore a real array, although it contains the
        !   information of a complex array.
        !
        ! Input
        ! =======
        ! - x : is the array to be  transformed
        ! - (m+2,n) : the shape of the array x;
        !               that is, the original array,
        !                padded with two rows.
        ! - plan : plan for the transform initialized by
        !               init_rfft2_ip
        ! Output
        ! ========
        ! - xhat : the complex-valued DFT of x

        implicit none

        integer, intent(in) :: m, n
        real (kind=8), dimension(m+2,n), intent(inout) :: x
        integer (kind=8), intent(in) :: plan

        ! compute forward 2d fft for real input
        call dfftw_execute_dft_r2c(plan,x,x)

    end subroutine rfft2_ip

    ! ifft2_2
    subroutine irfft2_ip(xhat,m,n,plan)

        ! It computes a 2D complex to real forward 2D fft
        !   in-place. The output overwrites the input.
        !
        ! Input
        ! =======
        ! - xhat : is the array to be  transformed (output of rfft2_ip)
        ! - (m+2,n) :    the shape of the array
        !
        ! - plan : plan for the transform initialized by
        !               init_irfft2_ip
        ! Output
        ! ========
        ! - x : the iDFT of xhat

        implicit none

        integer, intent(in) :: m,n
        real (kind=8), dimension(m+2,n), intent(inout) :: xhat
        integer (kind=8), intent(in) :: plan

        ! compute forward 1d fft for real input
        call dfftw_execute_dft_c2r(plan,xhat,xhat)
        xhat = xhat/n/n

    end subroutine irfft2_ip


    ! Wavenumber out-of-place
    subroutine rfft2_wavenumber(m,n,k,l)

        ! *** NOT WORKING YET, PLEASE BE CAREFUL... *** !


        ! It return the wavenumber arrays
        !  for rfft2 out-of-place

        implicit none

        integer, intent(in) :: m,n
        real (kind=8), dimension(m/2+1,n), intent(inout) :: k,l

        integer :: i, j

        do j=1,n/2
            k(:,j) = j-1.d0
            k(:,j+n/2) = j-n/2-1.d0
        end do

        do i=1,n/2+1
          l(i,:) = i-1.d0
        end do

    end subroutine rfft2_wavenumber

! Wavenumber in-of-place
subroutine rfft2_wavenumber_ip(m,n,k,l)

    ! It return the wavenumber arrays
    !  for rfft2 out-of-place

    implicit none

    integer, intent(in) :: m,n
    real (kind=8), dimension(m+2,n), intent(inout) :: k,l

    integer :: i, j

    do j=1,n/2
        k(:,j) = j-1
        k(:,j+n/2) = j-n/2-1
    end do

    do i=1,n/2+1
      l(i,:) = i-1
      l(i+n/2,:) = n/2-(i-1)
    end do

end subroutine rfft2_wavenumber_ip


end module fft_mod
