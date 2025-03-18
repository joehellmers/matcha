! Copyright (c), The Regents of the University of California
! Terms of use are as specified in LICENSE.txt
module matcha_mimetic_test_m
  use julienne_m, only : test_t, test_result_t
  use input_m, only : input_t
  use output_m, only : output_t
  use matcha_m, only : matcha
  implicit none
  
  private
  public :: matcha_mimetic_test_t
  
  type, extends(test_t) :: matcha_mimetic_test_t
  contains
    procedure, nopass :: subject
    procedure, nopass :: results
  end type
  
contains 

  pure function subject() result(specimen)
    character(len=:), allocatable :: specimen
    specimen = "A matcha_mimetic_t object"
  end function
  
  function results() result(test_results)
    type(test_result_t), allocatable :: test_results(:)
    
    character(len=*), parameter :: longest_description = &
      "creates a simulated distribution matching an empirical distribution on each image" 
   

    test_results = test_result_t( &
      [ character(len=len(longest_description)) :: &
        "mimetic operator testing" &
      ], &
      [ mimetic_testing() &
      ] &
    )
  end function
  
  function mimetic_testing() result(test_passes)
    logical test_passes

    ! Hard coding order (k) to 2 for now
    ! These will become parameters later

    double precision :: a = 0.0
    double precision :: b = 1.0

    integer :: k = 2        ! Order of accuracy
    integer :: m            ! Number of cells
    double precision :: dx        ! Cell spacing
    double precision :: d = 1.0   ! Dirichlet coefficient
    double precision :: n = 1.0   ! Neumann coefficient

    ! Indexes
    integer :: i

    !Mimetic operator allocation
    double precision, allocatable, dimension(:,:) :: mimeticDiv1D, mimeticDiv1DTest
    double precision, allocatable, dimension(:,:) :: mimeticGrad1D
    double precision, allocatable, dimension(:,:) :: mimeticLaplace1D

    ! Weights
    double precision, allocatable, dimension(:) :: Q
    double precision, allocatable, dimension(:) :: P

    ! Robin Boundary Conditions
    double precision, allocatable, dimension(:,:) :: A_RBC
    double precision, allocatable, dimension(:,:) :: BG_RBC
    double precision, allocatable, dimension(:,:) :: Grad1D_RBC
    double precision, allocatable, dimension(:,:) :: BC_RBC

    ! Staggered Grid
    double precision, allocatable, dimension(:) :: grid

    ! RHS
    double precision, allocatable, dimension(:) :: rhs

    test_passes = .true.

    m = 2*k + 1
    dx = (b-a)/m

    allocate(Q(m))
    allocate(P(m))

    Q = 1.0
    P = (/ (3.0/8.0), (9.0/8.0), 1.0, (9.0/8.0), (3.0/8.0) /)

    write(*,*) "Weights - Q:" 
    write(*,*) Q
    write(*,*) "Weights - P:" 
    write(*,*) P

    write(*,*) "Starting with 1-D Divergence"
    call makeMimeticDiv1D(mimeticDiv1D, k)
    write(*,25) mimeticDiv1D

    write(*,*) "Now do 1-D Gradient"
    call makeMimeticGrad1D(mimeticGrad1D, k)
    write(*,25) mimeticGrad1D

    write(*,*) "Now do 1-D Laplacian"
    call makeMimeticLaplace1D(mimeticLaplace1D, k)
    write(*,26) mimeticLaplace1D

    write(*,*) "Add 1-D Robin BCs"
    call addRobinBCs(mimeticLaplace1D, k, a, b)
    write(*,26) mimeticLaplace1D


    ! write(*,*) "Final Robin BC:"
    ! write(*,26) BC_RBC

    ! mimeticLaplace1D = mimeticLaplace1D + BC_RBC

    ! write(*,*) "Final Laplacian:"
    ! write(*,26) mimeticLaplace1D

    allocate(grid(m + 2))
    grid(1) = a;
    grid(2) = a + dx / 2.0
    do i = 3, m + 1
        grid(i) = grid(i - 1) + dx
    end do
    grid(m + 2) = b

    allocate(rhs(m + 2))
    rhs = exp(grid)
    rhs(1) = 0
    rhs(m + 2) = 2 * exp(1.0)

    deallocate(mimeticDiv1D)
    deallocate(mimeticGrad1D)
    deallocate(mimeticLaplace1D)
    ! deallocate(A_RBC)
    ! deallocate(BG_RBC)
    ! deallocate(BC_RBC)
    deallocate(grid)
    deallocate(rhs)


    ! Formatting
    25 format (6f7.2)
    26 format (7f7.2)

  end function

  subroutine makeMimeticDiv1D(DivOp1D, order)
    implicit none

    double precision, allocatable, intent(inout) :: DivOp1D(:,:)
    integer, intent(in) :: order

    integer :: m
    integer :: i
    
    m = order * 2 + 1

    allocate(DivOp1D(m+2,m+1))

    DivOp1D(1,1) = 0.0
    DivOp1D(m+2,m+1) = 0.0

    do i = 2, m+1
        DivOp1D(i,i-1) = -1.0
        DivOp1D(i,i) = 1.0;    
    end do

  end subroutine makeMimeticDiv1D


  subroutine makeMimeticGrad1D(GradOp1D, order)

    implicit none

    double precision, allocatable, intent(inout) :: GradOp1D(:,:)
    integer, intent(in) :: order

    integer :: m
    integer :: i
    
    m = order * 2 + 1
    allocate(GradOp1D(m+1,m+2))
    GradOp1D = 0.0

    GradOp1D(1,1) = -8.0/3.0
    GradOp1D(1,2) = 3.0
    GradOp1D(1,3) = -1.0/3.0

    GradOp1D(m+1,m+2) = 8.0/3.0
    GradOp1D(m+1,m+1) = -3.0
    GradOp1D(m+1,m) = 1.0/3.0

    do i = 2, m
        GradOp1D(i,i) = -1.0
        GradOp1D(i,i+1) = 1.0;    
    end do

  end subroutine makeMimeticGrad1D

  subroutine makeMimeticLaplace1D(LaplaceOp1D, order)

    implicit none

    double precision, allocatable, intent(inout) :: LaplaceOp1D(:,:)
    integer, intent(in) :: order

    double precision, allocatable :: DivOp1D(:,:)
    double precision, allocatable :: GradOp1D(:,:)

    integer :: m
    integer :: i
    
    m = order * 2 + 1

    allocate(LaplaceOp1D(m+2,m+2))

    call makeMimeticDiv1D(DivOp1D, order)
    call makeMimeticGrad1D(GradOp1D, order)

    LaplaceOp1D = matmul(DivOp1D,GradOp1D)

  end subroutine makeMimeticLaplace1D

  subroutine addRobinBCs(LaplacianOp1D, order, a, b)

    double precision, allocatable, intent(inout) :: LaplacianOp1D(:,:)
    integer, intent(in) :: order
    double precision, intent(in) :: a, b ! End points of 1-D domain

    double precision, allocatable :: A_RBC(:,:), BG_RBC(:,:), BC_RBC(:,:), Grad1D_RBC(:,:)

    integer :: m

    m = order * 2 + 1

    allocate(A_RBC(m+2,m+2))
    allocate(BG_RBC(m+2,m+2))
    allocate(BC_RBC(m+2,m+2))
    A_RBC = 0.0; BG_RBC = 0.0; BC_RBC = 0.0;


    A_RBC(1, 1) = a;
    A_RBC(m + 2, m + 2) = a;

    call makeMimeticGrad1D(Grad1D_RBC, order)  

    BG_RBC(1,:)         = -b * Grad1D_RBC(1,:);
    BG_RBC(m + 2, :)    =  b * Grad1D_RBC(m + 1, :);
    
    BC_RBC = A_RBC + BG_RBC

    LaplacianOp1D = LaplacianOp1D + BC_RBC

    deallocate(A_RBC)
    deallocate(BG_RBC)
    deallocate(BC_RBC)
    deallocate(Grad1D_RBC)

  end subroutine addRobinBCs

end module matcha_mimetic_test_m
