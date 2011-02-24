module line
  !>ANYTHING THAT ALTERS THE VORTEX FILAMENT (GEOMETRICALLY) SHOULD BE HERE
  use cdata
  use general
  use periodic
  contains
  subroutine pinsert
    !insert particles to maintain a ~constant resolution
    implicit none    
    type(qvort), allocatable, dimension(:) :: tmp
    real :: disti, curv, f_ddot(3)
    integer :: par_new
    integer :: old_pcount, i
    !check that we have particles
    if (quasip_only.eqv..false.) then
      if (count(mask=f(:)%infront>0)==0) then
        !call fatal_error('line.mod:pinsert','vortex line length is 0, run over')
      end if
    end if
    old_pcount=pcount
    total_length=0. !zero this
    if (magnetic) Brms=sqrt(sum(f(:)%B**2,mask=f(:)%infront>0)/count(mask=f(:)%infront>0)) !rms of magnetic field
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
        !average the current velocity
        f(par_new)%u=0.5*(f(i)%u+f(f(i)%infront)%u) 
        !zero the older velocities
        f(par_new)%u1=0. ; f(par_new)%u2=0.

        !set correct infront and behinds & ghostzones
        f(par_new)%behind=i ; f(par_new)%infront=f(i)%infront
        call get_ghost_p(par_new,f(par_new)%ghosti, f(par_new)%ghostb) !periodic.mod
        f(f(i)%infront)%behind=par_new
        call get_ghost_p(f(i)%infront,f(f(i)%infront)%ghosti, f(f(i)%infront)%ghostb) !periodic.mod         
        f(i)%infront=par_new
        call get_ghost_p(i,f(i)%ghosti, f(i)%ghostb) !periodic.mod         
        !finally address the magnetic issue
        if (magnetic) then
          if (f(i)%B<5.*Brms) then
              f(par_new)%B=2*f(i)%B
              f(i)%B=2*f(i)%B
          else
            f(par_new)%B=f(i)%B
          end if
        end if 
      end if
    end do
    !calculate average separation of particles
    avg_sep=total_length/old_pcount
    !print*, maxloc(f(:)%B), maxval(f(:)%B), f(maxloc(f(:)%B))%Bstretch
  end subroutine
  !*************************************************************************
  subroutine premove
    !maximum curvature cutoff based on an equilateral triangle with side
    !of length delta is \sqrt(3)/delta. This routine will remove particles
    !whose curvature is too high
    implicit none
    real :: distii
    integer :: infront, tinfront
    integer :: i
    logical :: do_remove 
    do i=1, pcount
      if (f(i)%infront==0) cycle !empty particle
      do_remove=.false.
      !get the distance between the particle and the one twice infront
      distii=distf(i,f(f(i)%infront)%infront)
      !check curvature if required
      if (phonon_emission) then
        if (curvature(f(i)%infront)>phonon_percent*(2./delta)) then
          if (distii<2.*delta) then!we don't remove if separation is too large
            do_remove=.true.
          end if
        end if
      end if  
      if (distii<delta) then
        do_remove=.true.
      end if
      if (do_remove) then
        !print to file the curvature of this particle
        infront=f(i)%infront ; tinfront=f(f(i)%infront)%infront
        open(unit=56,file='./data/removed_curv.log',position='append')
          write(56,*) curvature(infront)
        close(56)
        !remove the particle at infront
        f(tinfront)%behind=i ; f(i)%infront=tinfront
        call clear_particle(infront) !general.mod
        !check the size of the new loop
        call loop_killer(i)
        remove_count=remove_count+1
        !halve the magnetic field - due to contraction
        if (magnetic) then
          f(i)%B=f(i)%B/2.
        end if 
      end if
    end do
  end subroutine
  !******************************************************************
  !>find the closest particle to i using N^2 operation
  !>we need to do particles on the boundary as well (periodic)
  subroutine pclose
    implicit none
    integer :: i, j
    real :: dist
    do i=1, pcount
      f(i)%closestd=10. !arbitrarily high
      do j=1, pcount
        if ((i/=j).and.(f(i)%infront/=j).and.(f(i)%behind/=j)) then
        !the above line ensures we do not reconnect with particles infront/behind
          dist=distf(i,j)
          if (dist<f(i)%closestd) then
           f(i)%closest=j
           f(i)%closestd=dist
          end if
        end if
      end do
    end do
  end subroutine
  !******************************************************************
  subroutine precon
    !THE ROUTINE THAT RECONNECTS FILAMENTS WHICH BECOME CLOSE
    implicit none
    real :: distr, min_distr !reconnection distances
    real :: dot_val, tangent1(3), tangent2(3) !used to determine if filaments parallel
    integer :: pari, parb, parii, parbb, parji, parjb !particles infront/behind
    integer :: par_recon !the particle we reconnect with
    integer :: i, j !we must do a double loop over all particles N^2
    logical :: same_loop
    do i=1, pcount
      if (f(i)%infront==0) cycle !empty particle
      pari=f(i)%infront ; parb=f(i)%behind !find particle infront/behind
      parii=f(pari)%infront ; parbb=f(parb)%behind !find particle twice infront/behind
      !now we determine if we can reconnect
      if ((f(i)%closestd<delta/2.).and.(f(i)%closestd>epsilon(1.))) then
        j=f(i)%closest
        !another saftery check
        if (j==pari) cycle ; if (j==parb) cycle ; if (j==0) cycle
        !these two could have reconnected earlier in this case j will be empty
        if (f(j)%infront==0) cycle
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
          if (recon_info) then !more reconnection information
            call same_loop_test(i,j,same_loop) !line.mod
            if (same_loop) then
              self_rcount=self_rcount+1 
            else
              vv_rcount=vv_rcount+1
            end if
          end if
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
  !******************************************************************
  subroutine precon2
    !THE ROUTINE THAT RECONNECTS FILAMENTS WHICH BECOME CLOSE
    implicit none
    real :: distr, min_distr !reconnection distances
    real :: dot_val, tangent1(3), tangent2(3) !used to determine if filaments parallel
    integer :: pari, parb, parii, parbb, parji, parjb !particles infront/behind
    integer :: par_recon !the particle we reconnect with
    integer :: i, j !we must do a double loop over all particles N^2
    logical :: same_loop
    do i=1, pcount
      if (f(i)%infront==0) cycle !empty particle
      pari=f(i)%infront ; parb=f(i)%behind !find particle infront/behind
      parii=f(pari)%infront ; parbb=f(parb)%behind !find particle twice infront/behind
      !now we determine if we can reconnect
      if ((f(i)%closestd<delta/2.).and.(f(i)%closestd>epsilon(1.))) then
        j=f(i)%closest
        !another saftery check
        if (j==pari) cycle ; if (j==parb) cycle ; if (j==0) cycle
        !these two could have reconnected earlier in this case j will be empty
        if (f(j)%infront==0) cycle
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
          if (recon_info) then !more reconnection information
            call same_loop_test(i,j,same_loop) !line.mod
            if (same_loop) then
              self_rcount=self_rcount+1 
            else
              vv_rcount=vv_rcount+1
            end if
          end if
          !set correct behind_infront
          f(parjb)%infront=pari
          f(pari)%behind=parjb
          f(i)%infront=j
          f(j)%behind=i
          !check the size of these new loops
          call loop_killer(pari) ; call loop_killer(i)
        end if 
      end if
    end do
  end subroutine
  !**************************************************
  subroutine loop_killer(particle)
    !removes loops with less than three particles
    !this is needed to ensure derivatives can be calculated correctly
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
    if (counter<6) then
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
  !**************************************************
   subroutine same_loop_test(i,j,same_loop)
     use Cdata
     implicit none
     integer,intent(IN) :: i,j
     integer :: k
     integer :: next
     logical :: same_loop
     !aim  of routine is to find out wether the links are on the same loop
     same_loop=.false. !initial condition now try and prove if true
     next=i
     do k=1, pcount
       next=f(k)%infront
       if (next==j) then
         same_loop=.true.
         exit
       end if
       if (next==i) exit
     end do
   end subroutine
end module
