module line
  use cdata
  use general
  use periodic
  contains
  !**************************************************************
  subroutine pinsert
    implicit none    
    type(qvort), allocatable, dimension(:) :: tmp
    real :: disti, curv, f_ddot(3)
    integer :: par_new
    integer :: old_pcount, i
    old_pcount=pcount
    total_length=0. !zero this
    do i=1, old_pcount
      if (i>size(f)) then
        !bug check
        print*, 'I think there is a problem' ; exit
      end if
      if (f(i)%infront==0) cycle !empty particle
      !get the distance between the particle and the one infront
      disti=dist_gen(f(i)%x,f(i)%ghosti) !general.f90
      total_length=total_length+disti !measure total length of filaments
      if (disti>delta) then               
        !print*, 'in here' ; stop
        !we need a new particle 
        !the first step is assess where to put the particle?
        !1. is there an empty slot in out array?
        if (minval(f(:)%infront)==0) then
          !there is an empty slot in the array
          !now find its location minloc is fortran intrinsic function
          par_new=minloc(f(:)%infront,1)
        else
          !2. we must resize the array - 
          !increment by 1 (speed) - use a dummy array tmp
          allocate(tmp(size(f)+1)) ; tmp(:)%infront=0 !0 the infront array
          !copy accross information
          tmp(1:size(f)) = f
          !now deallocate tmp and transfer particles back to f
          !move_alloc is new intrinsic function (fortran 2003)
          call move_alloc(from=tmp,to=f)
          !pcount+1 must be free - so put the particle there
          par_new=pcount+1 ; pcount=pcount+1 !increase the particle count
        end if
        !insert a new particle between i and infront using curvature
        !get second derivative at i
        call get_deriv_2(i,f_ddot) !general.mod
        curv=sqrt(dot_product(f_ddot,f_ddot)) !get the curvature
        f(par_new)%x=0.5*(f(i)%x+f(i)%ghosti)+&
                     (sqrt(curv**2-0.25*disti**2)-curv)*curv*f_ddot
        !f(par_new)%x=(f(i)%x+f(i)%ghosti)/2.
        !0 the velocities
        f(par_new)%u=0. ; f(par_new)%u1=0. ; f(par_new)%u2=0.

        !set correct infront and behinds & ghostzones
        f(par_new)%behind=i ; f(par_new)%infront=f(i)%infront
        call get_ghost_p(par_new,f(par_new)%ghosti, f(par_new)%ghostb) !periodic.mod
        f(f(i)%infront)%behind=par_new
        call get_ghost_p(f(i)%infront,f(f(i)%infront)%ghosti, f(f(i)%infront)%ghostb) !periodic.mod         
        f(i)%infront=par_new
        call get_ghost_p(i,f(i)%ghosti, f(i)%ghostb) !periodic.mod         
      end if
    end do
    !calculate average separation of particles
    avg_sep=total_length/old_pcount
  end subroutine
  !******************************************************************
  subroutine premove
    !maximum curvature cutoff based on an equilateral triangle with side
    !of length delta is \sqrt(3)/delta. This routine will remove particles
    !whose curvature is too high
    implicit none
    real :: distii
    integer :: infront, tinfront
    integer :: i    
    do i=1, pcount
      if (f(i)%infront==0) cycle !empty particle
      !get the distance between the particle and the one twice infront
      distii=distf(i,f(f(i)%infront)%infront)
      if (distii<delta) then
        infront=f(i)%infront ; tinfront=f(f(i)%infront)%infront
        !remove the particle at infront
        f(tinfront)%behind=i ; f(i)%infront=tinfront
        call clear_particle(infront) !general.mod
      end if
      !check the size of the new loop
      call loop_killer(i)
    end do
  end subroutine
  !******************************************************************
  subroutine precon
    !THE ROUTINE THAT RECONNECTS FILAMENTS WHICH BECOME CLOSE
    !THIS IS $N^2$ OPERATION.
    implicit none
    real :: distr, min_distr !reconnection distances
    real :: dot_val, tangent1(3), tangent2(3) !used to determine if filaments parallel
    integer :: pari, parb, parii, parbb, parji, parjb !particles infront/behind
    integer :: par_recon !the particle we reconnect with
    integer :: i, j !we must do a double loop over all particles N^2
    do i=1, pcount
      if (f(i)%infront==0) cycle !empty particle
      min_distr=10. !arbitrarily large
      pari=f(i)%infront ; parb=f(i)%behind !find particle infront/behind
      parii=f(pari)%infront ; parbb=f(parb)%behind !find particle twice infront/behind
      do j=i, pcount !start the seocnd loop (running upwards from i is sufficient)
        if ((i/=j).and.(pari/=j).and.(parb/=j).and.(parii/=j).and.(parbb/=j)) then
          !the above line ensures we do not reconnect with particles infront/behind
          !and also particles twice infront/behind, these operations are accounted for
          !in pinsert/premove

          !get the distance between i and j
          distr=distf(i,j) !general.mod
          if (distr<min_distr) then
            !store this distance/particle
            min_distr=distr ; par_recon=j
          end if
        end if
      end do
      !now we determine if we can reconnect
      if (min_distr<delta/2.) then
        j=par_recon !go back to working with j (shorter!)
        parji=f(j)%infront ; parjb=f(j)%behind
        !we can reconnect based on distance
        !now check whether parallel
        tangent1=norm_tanf(i) ; tangent2=norm_tanf(j) !general.mod
        dot_val=dot_product(tangent1,tangent2) !intrinsic function
        if ((dot_val>0.9)) then
          !we cannot reconnect as filaments are parallel          
          !print*, 'cannot reconnect, cos(theta) is:',dot_val
        else
          !print*, 'i',i, f(i)%infront, f(i)%behind
          !print*, 'j',j, f(j)%infront, f(j)%behind
          !reconnect the filaments
          recon_count=recon_count+1 !keep track of the total # of recons
          !reomove two particles involved in reconnection
          call clear_particle(i) ; call clear_particle(j)
          !set correct behind_infront
          f(parjb)%infront=pari
          f(pari)%behind=parjb
          f(parb)%infront=parji
          f(parji)%behind=parb
          !check the size of these new loops
          call loop_killer(pari) ; call loop_killer(parb)
        end if 
      end if
    end do
  end subroutine
  !**************************************************
  subroutine loop_killer(particle)
    !removes loops with less than three particles
    implicit none
    integer :: particle, next
    integer :: store_next
    integer :: i, counter
    counter=1
    next=particle
    do i=1, pcount            
      next=f(next)%infront
      if (next==particle) exit  
      counter=counter+1
    end do
    ! If loop is too small destroy
    if (counter<4) then
      next=particle 
      do i=1, pcount
        store_next=f(next)%infront
        call clear_particle(next) !general.mod
        next=store_next
        if (next==particle) then
          exit  
        end if
      end do
    end if
  end subroutine
end module
