! file: fast_logdet.f90
! compile with f2py -c -m fortran_v fortran_version.f90

subroutine vec_logdet(x, m, n, logdets)
    implicit none
    integer m, n, i, j, info
    !f2py integer intent(hide),depend(x) :: n = shape(x,0), m=shape(x,2)
    !m and n should be equal
    !make sure to pass in X as column major for convenience
    real*8, intent(in):: x(n,n,m)
    real*8, intent(out):: logdets(m)
    do i=1, m
        logdets(i) = 0d0
        call dpotrf('U', n, x(:,:,i), n, info)
        do j=1, n
             logdets(i) = logdets(i) + log(x(j,j,i))
        end do
        logdets(i) = 2d0*logdets(i)
    end do
end

subroutine fast_logdet(x, n, logdet)
    implicit none
    integer n, j, info
    !f2py integer intent(hide),depend(x) :: n = shape(x,0)
    !make sure to pass in X as column major for convenience
    !also note that this operation is in place
    real*8, intent(in):: x(n,n)
    real*8, intent(out):: logdet
    logdet = 0d0
    call dpotrf('U', n, x, n, info)
    if (info /= 0) then
        print *, "matrix not positive definite"
    end if
    do j=1, n
        logdet = logdet + log(x(j,j))
    end do
    logdet = 2d0*logdet
end

!subroutine best_from_many_outer_rot_roundings(Sigma, X, n_Gammas, k)
!    implicit none
!end 
!
!subroutine many_outer_rot_roundings(Sigma, X, n_Gammas)
!    implicit none
!end
!
!subroutine looped_outer_rotation(Gammas, U)
!    implicit none
!end

subroutine find_best_in_subtree(Sigma, u, n, r, nz, u_best, best_obj, leaf_nodes)
    !the 'all_combinations' part of this code comes from 
    !https://docs.python.org/3/library/itertools.html#itertools.combinations_with_replacement
    !k is the toal number of sensors to place.
    !r is the number of sensors left to place.
    !u is the placements
    implicit none
    integer n, ix, i,j
    !f2py integer intent(hide),depend(u) :: n = shape(u,0)
    integer, intent(in):: r !number of sensors left to place
    integer, intent(in):: nz !number of undecided candidate locations
    integer, intent(out):: leaf_nodes !number of leaf nodes visited
    real*8, intent(in):: u(n)
    real*8, intent(in):: Sigma(n,n)
    real*8 Sigma_cp(n,n)
    real*8 u_cp(n)
    real*8, intent(out):: u_best(n)
    real*8, intent(out):: best_obj
    real*8 cur_obj, logdet_sigma
    integer zero_inds(nz)
    integer selected_inds(r)
    j=1 !counter for nonzero indices
    do i=1,n
        if (u(i)==0d0) then
            zero_inds(j) = i
            j = j+1
        end if
    end do
    selected_inds = [(ix, ix=1,r)]
    u_cp = u
    u_cp(zero_inds) = -1d0
    !call choose(n, k, num_comb) 
    !do first calculation here. Then we move on to the next calculation
    u_cp(zero_inds(selected_inds)) = 1d0
    !print *,zero_inds
    !print *,zero_inds(selected_inds)
    !print *,u_cp
    call obj(Sigma, u_cp, n, best_obj)
    leaf_nodes = leaf_nodes + 1
    u_best = u_cp !start with this value and then update 
    do while (0 < 1)
        ix = 0
        do i=r, 1, -1
            if (selected_inds(i) /= i+nz-r) then !TODO: Check this logic
                ix = 1
                exit
            end if
        end do
        if (ix == 0) then !we've tried all placements
            Sigma_cp = Sigma
            !copy because this operation is in place
            call fast_logdet(Sigma_cp, n, logdet_sigma)
            best_obj = best_obj - logdet_sigma
            return
        end if
        selected_inds(i) = selected_inds(i) + 1
        do j=i+1,r
            selected_inds(j) = selected_inds(j-1) + 1 
        end do        
        u_cp(zero_inds) = -1d0
        u_cp(zero_inds(selected_inds)) = 1d0
        !print *,zero_inds
        !print *,zero_inds(selected_inds)
        !print *,u_cp
        call obj(Sigma, u_cp, n, cur_obj)
        leaf_nodes = leaf_nodes + 1
        if (cur_obj > best_obj) then
            best_obj = cur_obj
            u_best = u_cp
        end if
   end do
end

subroutine obj(Sigma, u, n, logdet)
    implicit none
    integer i,j
    integer, intent(in):: n
    real*8, intent(in):: u(n)
    real*8, intent(in):: Sigma(n,n)
    real*8, intent(out):: logdet
    real*8 Z(n,n)
    do j=1,n
        do i=1,n
            Z(i,j) = Sigma(i,j)*(u(i)*u(j) + 1d0)/2d0
        end do
    end do
    call fast_logdet(Z, n, logdet) 
end 

subroutine spm_greedy(Sigma, n, k, greedy_u, greedy_obj)
    implicit none
    integer i,j
    integer, intent(in):: n
    !f2py integer intent(hide),depend(Sigma) :: n = shape(Sigma,0)
    integer, intent(in):: k
    real*8, intent(in):: Sigma(n,n)
    real*8 Sigma_cp(n,n)
    real*8, intent(out):: greedy_u(n)
    real*8, intent(out):: greedy_obj
    real*8 u_cp(n)
    real*8 cur_obj
    integer max_ind
    max_ind = 1 !avoids a compiler warning
    greedy_u = -1d0
    do i=1,k
        do j=1,n
            u_cp = greedy_u
            if (u_cp(j) /= 1d0) then
               u_cp(j) = 1d0
               call obj(Sigma, u_cp, n, cur_obj)
               if (j == 1 .or. cur_obj > greedy_obj) then
                  greedy_obj = cur_obj
                  max_ind = j
               end if
            end if
        end do
        greedy_u(max_ind) = 1d0
    end do
    Sigma_cp = Sigma
    !copy because this operation is in place
    !we just reuse cur_obj variable as the logdet of Sigma
    call fast_logdet(Sigma_cp, n, cur_obj)
    greedy_obj = greedy_obj - cur_obj
end

subroutine choosebranchingvar_greedy_good(Sigma, n, u, ind)
    implicit none
    integer j, num_unplaced
    integer, intent(in):: n
    !f2py integer intent(hide),depend(u) :: n = shape(u,0)
    real*8, intent(in):: Sigma(n,n)
    real*8, intent(in):: u(n)
    real*8 u_cp(n)
    real*8 cur_obj
    real*8 best_obj
    integer, intent(out):: ind
    best_obj = -1e6 !get started with a large negative number
    u_cp = u
    do j=1,n
        if (u_cp(j)==0d0) then
            u_cp(j) = -1d0
        end if
    end do
    num_unplaced = 0
    do j=1,n
        if (u(j) == 0d0) then
           u_cp(j) = 1d0
           call obj(Sigma, u_cp, n, cur_obj)
           if (num_unplaced==0 .or. cur_obj > best_obj) then !j == 1 .or. 
              best_obj = cur_obj ! .or. cur_obj > best_obj
              ind = j
           end if
           u_cp(j) = -1d0 ! set it back to -1d0
           num_unplaced = num_unplaced + 1
        end if
    end do
    ind = ind - 1 !Python vs Fortran indexing
end

! bad greedy method, chooses worst objective
subroutine choosebranchingvar_greedy_bad(Sigma, n, u, ind)
    implicit none
    integer j, num_unplaced
    integer, intent(in):: n
    !f2py integer intent(hide),depend(u) :: n = shape(u,0)
    real*8, intent(in):: Sigma(n,n)
    real*8, intent(in):: u(n)
    real*8 u_cp(n)
    real*8 cur_obj
    real*8 best_obj
    integer, intent(out):: ind
    best_obj = 1e6 !get started with a large negative number
    u_cp = u
    do j=1,n
        if (u_cp(j)==0d0) then
            u_cp(j) = -1d0
        end if
    end do
    num_unplaced = 0
    do j=1,n
        if (u(j) == 0d0) then
           u_cp(j) = 1d0
           call obj(Sigma, u_cp, n, cur_obj)
           if (num_unplaced==0 .or. cur_obj < best_obj) then !j == 1 .or. 
              best_obj = cur_obj ! .or. cur_obj > best_obj
              ind = j
           end if
           u_cp(j) = -1d0 ! set it back to -1d0
           num_unplaced = num_unplaced + 1
        end if
    end do
    ind = ind - 1 !Python vs Fortran indexing
end

! choose the branching variable randomly
subroutine choosebranchingvar_random(n, u, ind)
    implicit none
    integer, intent(in):: n
    !f2py integer intent(hide),depend(u) :: n = shape(u,0)
    
    real*8, intent(in):: u(n)
    ! real*8 u_cp(n)
    ! integer max_ind
    ! real*8 max_Sigma
    integer, intent(out):: ind
    real*8 l

    !set index to a random number
    call random_number(l)
    ind = 1+FLOOR(n*l) 
    
    do while (u(ind) .ne. 0d0)  !if we've chosen an index which has already been placed, resample
        call random_number(l)
        ind = 1+FLOOR(n*l)
    end do

    ind = ind - 1 !Python vs Fortran indexing
end
    
!altered to choose max variance
subroutine choosebranchingvar_max_var(Sigma, n, u, ind)
    implicit none
    integer j
    integer, intent(in):: n
    !f2py integer intent(hide),depend(u) :: n = shape(u,0)
    real*8, intent(in):: Sigma(n,n)
    real*8, intent(in):: u(n)
    ! real*8 u_cp(n)
    ! integer max_ind
    real*8 max_Sigma
    integer, intent(out):: ind

    !set max varinace index as index 1
    !set to 0 for catching case where all are placed
    ind = 0
    !set max variance as Sigma(1,1)
    max_Sigma = 0

    !loop and check if any variance is larger
    do j=1,n
        
        if(u(j) == 0 .and. Sigma(j,j) > max_Sigma) then
            max_Sigma = Sigma(j,j)
            ind = j
        end if

    end do

    ind = ind - 1 
end

subroutine best_from_many_outer_rot_roundings(Sigma, n, X, n_Gammas, k, best_u, best_obj)
    implicit none
    integer, intent(in):: n
    !f2py integer intent(hide),depend(Sigma) :: n = shape(Sigma,0)
    real*8, intent(in):: Sigma(n,n)
    real*8 Sigma_cp(n,n)
    real*8, intent(in):: X(n,n)
    real*8 X_cp(n,n), eig_vals(n), eig_vecs(n,n), U(n,n)
    integer, intent(in):: n_Gammas
    integer, intent(in):: k
    integer num_ones, i, info, j, m, isuppz(2*n)
    real*8, intent(out)::best_obj
    real*8, intent(out)::best_u(n)
    real*8, parameter :: r8_pi = 3.141592653589793D+00
    real*8, allocatable :: work(:)
    integer, allocatable :: iwork(:)
    integer opt_iwork_size(1)
    real*8 opt_work_size(1)
    real*8 cur_u(n), p1(n), p2(n)
    real*8 cur_obj
    real*8 gam
    real*8 tol 
    tol = -1 
    X_cp = X
    best_obj = -1e6 !a placeholder in case there's nothing feasible
    !we call dpstrf instead of dpotrf b/c X_cp can be degenerate
    !note this operates in place, so we can 
    !note we take abstol=1e-6, the numerical
    !tolerance for finding eigenvalues.
    !could do something inspired by scs's
    !tolerance instead
    !first we query to see the optimal size of work and iwork
    !print *, "prior to finding optimal array sizes" 
    call dsyevr('V', 'A', 'L', n, X_cp, n, 0d0, 0d0, 0, 0, -1e6, m, & 
                eig_vals, eig_vecs, n, isuppz, opt_work_size, -1, &
                opt_iwork_size, -1, info)
    !allocate work and iwork according to their optimal sizes
    allocate(work(int(opt_work_size(1))), stat=info)
    allocate(iwork(opt_iwork_size(1)), stat=info)
    !do the eigendecomposition
    !print *, "optimal work size is", int(opt_work_size(1))
    !print *, "optimal iwork size is", opt_iwork_size(1)
    call dsyevr('V', 'A', 'L', n, X_cp, n, 0d0, 0d0, 0, 0, -1e6, m, &
                eig_vals, eig_vecs, n, isuppz, work, int(opt_work_size(1)), &
                iwork, opt_iwork_size(1), info)
    if (info /= 0) then
        print *, "There was an error in the eigendecomposition."
        print *, "info = ", info
    end if
    !threshold the eig_vals to be nonnegative and compute U
    do i=1,n
        if (eig_vals(i) < 0d0) then
            eig_vals(i) = 0d0
        end if
        do j=1,n
            U(i,j) = sqrt(eig_vals(i))*eig_vecs(j,i)
        end do
    end do

    do j=0,n_Gammas-1
        gam = r8_pi/2d0*dble(j)/dble(n_Gammas)
        call std_normal_vector(n, p1)
        call std_normal_vector(n, p2)
        !we do the rotation and rounding on the same piece of memory, cur_u
        call dgemv('T', n, n, 1d0, U, n, p1, 1, 0d0, p1, 1)
        cur_u = cos(gam)*p1 + sin(gam)*p2
        num_ones = 0
        do i=1,n
            if (cur_u(i) < 0d0) then
                cur_u(i) = 1d0
                num_ones = num_ones + 1
            else
                cur_u(i) = -1d0
            end if
        end do
        call obj(Sigma, cur_u, n, cur_obj)
        if (cur_obj > best_obj .and. (num_ones==k .or. num_ones==n-k)) then
            if (num_ones==n-k) then
                best_u = -1d0*cur_u
                best_obj = cur_obj
            else
                best_u = cur_u
                best_obj = cur_obj
            end if
         end if
    end do
    Sigma_cp = Sigma
    !copy because this operation is in place
    !we just reuse cur_obj variable as the logdet of Sigma
    call fast_logdet(Sigma_cp, n, cur_obj)
    best_obj = best_obj - cur_obj

    deallocate(work)
    deallocate(iwork)
end 

subroutine std_normal_vector(n, x)
  !modified from some code of John Burkardt, which is released
  !under LGPL license
  !https://people.sc.fsu.edu/~jburkardt/f_src/normal/normal.f90
  !Question: what to do about the seed?
  implicit none
  integer i
  integer, intent(in):: n !length of vector to return
  real*8 r1, r2
  real*8, parameter :: r8_pi = 3.141592653589793D+00
  real*8, intent(out):: x(n)
  do i=1,n
      call random_number ( harvest = r1 )
      call random_number ( harvest = r2 )
      x(i) = sqrt ( - 2.0D+00 * log ( r1 ) ) * cos ( 2.0D+00 * r8_pi * r2 )
  end do
end

subroutine choose(n, k, c)
!    ! taken from here: https://programming-idioms.org/idiom/67/binomial-coefficient-n-choose-k/3595/fortran
    implicit none
    integer, parameter :: i8 = selected_int_kind(18)
    integer, parameter :: dp = selected_real_kind(15)
!    !assumes k <= n
    integer, intent(in):: n, k
    integer*8, intent(out):: c
    c = nint(exp(log_gamma(n+1.0_dp)-log_gamma(n-k+1.0_dp)-log_gamma(k+1.0_dp)),kind=i8)
end



