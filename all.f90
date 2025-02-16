module distribution_m
  use iso_c_binding, only : c_double, c_int
  implicit none
  
  type distribution_t
    double precision, allocatable, dimension(:) :: vel_, cumulative_distribution_
  end type  

contains
  
  pure function construct(sample_distribution) result(distribution)
    double precision, intent(in) :: sample_distribution(:,:)
    type(distribution_t) distribution
    
    integer i
    if (.not. all(sample_distribution(:,2)>=0.D0)) error stop "negative sample_distribution value(s)"
    associate(nintervals => size(sample_distribution,1))      
      distribution%vel_ = [(sample_distribution(i,1), i =1, nintervals)]  ! Assign speeds to each distribution bin         
      distribution%cumulative_distribution_ = [0.D0, [(sum(sample_distribution(1:i,2)), i=1, nintervals)]]
      associate(f => distribution%cumulative_distribution_)
        if (.not. all([(f(i+1) >= f(i), i=1, size(f)-1)])) error stop "non-monotonic cum dist"
      end associate
    end associate
  end function

  pure function velocities(self, speeds, directions) result(my_velocities)
    class(distribution_t), intent(in) :: self
    double precision, intent(in) :: speeds(:,:), directions(:,:,:)
    double precision, allocatable :: my_velocities(:,:,:), sampled_speeds(:,:),  dir(:,:,:)
    
    if (.not. allocated(self%cumulative_distribution_)) error stop "unallocatd cum dist"
    if (.not. allocated(self%vel_)) error stop "unallocated vel_"

    call do_concurrent_sampled_speeds(speeds, self%vel_, self%cumulative_distribution_, sampled_speeds)

    associate(nsteps => size(speeds,2))
      dir = directions(:,1:nsteps,:)
      associate(dir_mag => sqrt(dir(:,:,1)**2 +dir(:,:,2)**2 + dir(:,:,3)**2))
        associate(dir_mag_ => merge(dir_mag, epsilon(dir_mag), dir_mag/=0.))
          dir(:,:,1) = dir(:,:,1)/dir_mag_
          dir(:,:,2) = dir(:,:,2)/dir_mag_
          dir(:,:,3) = dir(:,:,3)/dir_mag_
        end associate
      end associate
      call do_concurrent_my_velocities(nsteps, dir, sampled_speeds, my_velocities)
    end associate
  end function

  pure subroutine do_concurrent_sampled_speeds(speeds, vel, cumulative_distribution, sampled_speeds) bind(C)
    real(c_double), intent(in) :: speeds(:,:), vel(:), cumulative_distribution(:)
    real(c_double), intent(out), allocatable :: sampled_speeds(:,:)
    integer cell, step
    associate(ncells => size(speeds,1), nsteps => size(speeds,2))
      allocate(sampled_speeds(ncells,nsteps))
      do concurrent(cell = 1:ncells, step = 1:nsteps)
        associate(k => findloc(speeds(cell,step) >= cumulative_distribution, value=.false., dim=1)-1)
          sampled_speeds(cell,step) = vel(k)
        end associate
      end do
    end associate
  end subroutine
  
  pure subroutine do_concurrent_my_velocities(nsteps, dir, sampled_speeds, my_velocities) bind(C)
    integer(c_int), intent(in) :: nsteps
    real(c_double), intent(in) :: dir(:,:,:), sampled_speeds(:,:)
    real(c_double), intent(out), allocatable :: my_velocities(:,:,:)
    integer step
    
    if(allocated(my_velocities)) deallocate(my_velocities)
    allocate(my_velocities, mold=dir)
    
    do concurrent(step=1:nsteps)
      my_velocities(:,step,1) = sampled_speeds(:,step)*dir(:,step,1)
      my_velocities(:,step,2) = sampled_speeds(:,step)*dir(:,step,2)
      my_velocities(:,step,3) = sampled_speeds(:,step)*dir(:,step,3)
    end do
  end subroutine

end module distribution_m

  use distribution_m, only : distribution_t, velocities, construct
  implicit none
  integer, parameter  :: ncells = 6000, npositions = 6000, ndim = 3, nveldim = 4, nsteps = npositions - 1
  double precision random_4vectors(ncells,nsteps,nveldim)

  call random_init(repeatable=.true., image_distinct=.true.)
  call random_number(random_4vectors)  
  associate(v => velocities(construct(sample_distribution()), random_4vectors(:,:,1), random_4vectors(:,:,2:4)))
  end associate

contains

  function sample_distribution()
     double precision, allocatable :: sample_distribution(:,:), speeds(:), probability(:)
     double precision, parameter :: two_pi = 2D0*acos(-1.d0), speed_lower = 0.d0, speed_upper = 6.d0
     integer, parameter :: nintervals = 4
     integer i

     allocate(speeds(nintervals), probability(nintervals), sample_distribution(nintervals,2))
     
     associate(range => speed_upper - speed_lower)
       associate(dspeed => range/dble(nintervals))
        do i = 1,nintervals
          associate(speed_lower_bin => speed_lower + dble(i-1)*dspeed, speed_upper_bin => speed_lower + dble(i)*dspeed)
            speeds(i) = 0.5D0*(speed_lower_bin + speed_upper_bin)
          end associate
          probability(i) = exp(-(speeds(i)-3.d0)**2/2.d0)/dsqrt(two_pi) ! Use normal distribution
        end do
       end associate
     end associate

     sample_distribution(:,1) = speeds
     sample_distribution(:,2) = probability/sum(probability)
  end function

end
