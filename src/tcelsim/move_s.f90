submodule(move_m) t_cells
  implicit none

contains
  module procedure move_tcells
      ! Local variables
      integer, intent(out) i,j,k
      double precision, intent(out) speed,dt
      double precision, intent(out) rr1,rr2,rr3,sum

      ! Time step
      dt = .1

      do i = 1,ncells
        do j = 2,npositions

          ! Sample from the distribution
          call random_number(rr1)
          do k = 1,nintervals
            if (rr1 .ge. cumulative_distribution(k-1) .and. &
                rr1 .lt. cumulative_distribution(k)) then
               speed = vel(k)
            end if
          end do

          ! Create a random unit vector
          call random_number(rr1)
          call random_number(rr2)
          call random_number(rr3)
          sum = rr1 + rr2 + rr3
          rr1 = rr1/sum
          rr2 = rr2/sum
          rr3 = rr3/sum

          !     Use a forward Euler to advance the cell position
          x(i,j) = x(i,j-1) + dt*speed
          y(i,j) = y(i,j-1) + dt*speed
          z(i,j) = z(i,j-1) + dt*speed
        end do
      end do

    end module procedure move_tcells

end submodule move_m
