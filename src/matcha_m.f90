module matcha_m 
  use t_cell_collection_m, only : t_cell_collection_t
  use distribution_m, only : distribution_t
  use input_m, only : input_t
  implicit none

  interface matcha

    module function matcha(input) result(history)
      implicit none
      type(input_t), intent(in) :: input
      type(t_cell_collection_t), allocatable :: history(:)
    end function

  end interface
  
contains

  module procedure matcha
    double precision, parameter :: scale = 100.D0
    double precision, allocatable :: random_positions(:,:)

    associate( &
      ncells => input%num_cells(), &
      npositions => input%num_positions(), &
      ndim => input%num_dimensions(), &
      nintervals => input%num_intervals(), &
      dt => input%time_step() &
    )
    associate(nsteps => npositions -1)

    allocate(random_positions(ncells,ndim))
    call random_number(random_positions)  
    history = [t_cell_collection_t(scale*random_positions, time=0.D0)]

    block
      integer, parameter :: nveldim = 4
      double precision, allocatable :: random_4vectors(:,:,:), sample_distribution(:), velocities(:,:,:)
      type(distribution_t) distribution
      integer step
      
      allocate(random_4vectors(ncells,nsteps,nveldim))
      call random_number(random_4vectors)  
      allocate(sample_distribution(nintervals))
      call random_number(sample_distribution)
      sample_distribution = sample_distribution/sum(sample_distribution)
      distribution = distribution_t(sample_distribution)
      
      associate(random_speeds => random_4vectors(:,:,1), random_directions => random_4vectors(:,:,2:4))
        associate(v => distribution%velocities(random_speeds, random_directions))
          do step = 1, nsteps
            associate(x => history(step)%positions(), t => history(step)%time())
              history = [history, t_cell_collection_t(x + v(:,step,:)*dt, t + dt)]
            end associate
          end do
        end associate
      end associate
    end block

    end associate
    end associate

  end procedure

end module
