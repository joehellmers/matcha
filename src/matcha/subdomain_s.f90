submodule(subdomain_m) subdomain_s
  use data_partition_m, only : data_partition_t
  use assert_m, only : assert
  use intrinsic_array_m, only : intrinsic_array_t
  implicit none

  real, allocatable :: halo_x(:,:)[:]
  integer, parameter :: west=1, east=2

  type(data_partition_t) data_partition

contains

    module procedure define

      associate(nx => (n), ny => (n), num_subdomains => num_images(), me => this_image())

        call assert(num_subdomains <= nx, "subdomain_t%define: num_subdomains <= nx", intrinsic_array_t([nx, num_subdomains]))

        allocate(halo_x(ny,west:east)[*])
        call data_partition%define_partitions(nx)

        associate(my_first => data_partition%first(me), my_last => data_partition%last(me))
          allocate(self%s_(my_first:my_last, ny))
          block
            integer i

            self%s_(my_first:my_last, 1) = boundary_val
            do concurrent(i=2:ny-1)
              self%s_(my_first, i) = merge(boundary_val, internal_val, my_first==1)
              self%s_(my_first+1:my_last-1, i) = internal_val
              self%s_(my_last, i) = merge(boundary_val, internal_val, my_last==nx)
            end do
            self%s_(my_first:my_last, ny) = boundary_val
          end block
          self%dx_ = side/num_subdomains
          self%dy_ = self%dx_
        end associate
      end associate
    end procedure

    module procedure dx
      my_dx = self%dx_
    end procedure

    module procedure dy
      my_dy = self%dy_
    end procedure

    module procedure laplacian
      allocate(laplacian_rhs%s_, mold=rhs%s_)
      laplacian_rhs%s_ = 0.

      !laplacian_s%f = ( &
      !    s%f(1:nx-2,2:ny-1) - 2*s%f(2:nx-1, 2:ny-1) + f(3:nx, 2:ny-1) + &
      !    s%f(2:nx-1,1:ny-2) - 2*s%f(2:nx-1, 2:ny-1) + f(2:nx-1, 3:ny) &
      ! ) / dy**2
    end procedure

    module procedure multiply
      product%s_ =  lhs * rhs%s_
    end procedure

    module procedure add
      total%s_ =  lhs%s_ + rhs%s_
    end procedure

end submodule subdomain_s
