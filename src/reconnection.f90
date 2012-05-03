!>all reconnection routines are contained in this module this then feeds into line.mod
module reconnection
  use cdata
  use general
  implicit none
  contains
  !******************************************************************
  !>reconnect filaments if they become too close removing the two points
  !>which are reconnected, the two closest points. This is a very dissipative 
  !>scheme
  subroutine precon_dissapitive
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
        !we can reconnect based on distancerecon
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
            call same_loop_test(i,j,same_loop) !general.mod
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
          !set reconnection times
          f(parjb)%t_recon(2)=f(parjb)%t_recon(1) !move last recon to slot 2
          f(parjb)%t_recon(1)=t !set current time
          f(pari)%t_recon(2)=f(pari)%t_recon(1) 
          f(pari)%t_recon(1)=t 
          f(parb)%t_recon(2)=f(parb)%t_recon(1) 
          f(parb)%t_recon(1)=t 
          f(parji)%t_recon(2)=f(parji)%t_recon(1) 
          f(parji)%t_recon(1)=t 
          !check the size of these new loops
          call loop_killer(pari) ; call loop_killer(parb)
        end if 
      end if
    end do
  end subroutine
  !******************************************************************
  !>reconnect filaments if they become too close removing the two points
  !>which are reconnected, the two closest points. This is the scheme of
  !>Schwarz
  subroutine precon_schwarz
    implicit none
    real :: distr, min_distr !reconnection distances
    real :: dot_val, tangent1(3), tangent2(3) !used to determine if filaments parallel
    real :: schwarz_distance, mean_curv
    integer :: pari, parb, parii, parbb, parji, parjb !particles infront/behind
    integer :: par_recon !the particle we reconnect with
    integer :: i, j !we must do a double loop over all particles N^2
    logical :: same_loop
    do i=1, pcount
      if (f(i)%infront==0) cycle !empty particle
      pari=f(i)%infront ; parb=f(i)%behind !find particle infront/behind
      parii=f(pari)%infront ; parbb=f(parb)%behind !find particle twice infront/behind
      mean_curv=1./curvature(i)
      if (mean_curv>epsilon(0.)) then
        schwarz_distance=2.*mean_curv/(log(mean_curv/corea))
      else
        schwarz_distance=0.
      end if
      !now we determine if we can reconnect
      if ((f(i)%closestd<schwarz_distance).and.(f(i)%closestd>epsilon(1.))) then
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
        else
          !reconnect the filaments
          recon_count=recon_count+1 !keep track of the total # of recons
          if (recon_info) then !more reconnection information
            call same_loop_test(i,j,same_loop) !general.mod
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
          !set reconnection times
          f(parjb)%t_recon(2)=f(parjb)%t_recon(1) !move last recon to slot 2
          f(parjb)%t_recon(1)=t !set current time
          f(pari)%t_recon(2)=f(pari)%t_recon(1) 
          f(pari)%t_recon(1)=t 
          f(parb)%t_recon(2)=f(parb)%t_recon(1) 
          f(parb)%t_recon(1)=t 
          f(parji)%t_recon(2)=f(parji)%t_recon(1) 
          f(parji)%t_recon(1)=t 
          !check the size of these new loops
          call loop_killer(pari) ; call loop_killer(parb)
        end if 
      end if
    end do
  end subroutine
  !******************************************************************
  !>reconnect filaments if they become too close - default scheme
  subroutine precon_original
    implicit none
    real :: distr, min_distr !reconnection distances
    real :: dot_val, tangent1(3), tangent2(3) !used to determine if filaments parallel
    real :: l_before, l_after !to check line length before and after
    integer :: arc_dist !separation of points along the filament
    integer :: pari, parb, parii, parbb, parji, parjb !particles infront/behind
    integer :: par_recon !the particle we reconnect with
    integer :: i, j, k !we must do a double loop over all particles N^2
    logical :: same_loop
    do i=1, pcount
      if (f(i)%infront==0) cycle !empty particle
      do k=1, n_recon_track
        if (full_recon_distance(k)%active) then
          !do not reconnect points being monitored for reconnection dist
          !this is probably slightly dodgy!
          if (i==full_recon_distance(k)%i) cycle
          if (i==full_recon_distance(k)%j) cycle
        end if
      end do
      pari=f(i)%infront ; parb=f(i)%behind !find particle infront/behind
      parii=f(pari)%infront ; parbb=f(parb)%behind !find particle twice infront/behind
      !now we determine if we can reconnect
      if ((f(i)%closestd<delta/2.).and.(f(i)%closestd>epsilon(1.))) then
        j=f(i)%closest
        do k=1, n_recon_track
          if (full_recon_distance(k)%active) then
            !do not reconnect points being monitored for reconnection dist
            !this is probably slightly dodgy!
            if (j==full_recon_distance(k)%i) cycle
            if (j==full_recon_distance(k)%j) cycle
          end if
        end do
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
        else
          !now we must check the line length before and after
          l_before=dist_gen(f(i)%x,f(i)%ghosti)+dist_gen(f(i)%x,f(i)%ghostb)+&
                   dist_gen(f(j)%x,f(j)%ghosti)+dist_gen(f(j)%x,f(j)%ghostb)
          l_after=dist_gen(f(i)%x,f(i)%ghostb)+dist_gen(f(i)%x,f(j)%x)+&
                   dist_gen(f(j)%x,f(j)%ghosti)+dist_gen(f(pari)%x,f(parjb)%x)
          if (l_after<=l_before) then
            !reconnect the filaments
            recon_count=recon_count+1 !keep track of the total # of recons
            if (recon_info) then !more reconnection information
              call same_loop_test_arc_length(i,j,same_loop,arc_dist) !general.mod
              if (same_loop) then
                self_rcount=self_rcount+1 
                self_r_length=self_r_length+(l_before-l_after)
                open(unit=61,file='./data/recon_same_loop_sep.log',position='append')
                  write(61,*) t, arc_dist*avg_sep, arc_dist
                close(61)
              else
                vv_rcount=vv_rcount+1
                vv_r_length=vv_r_length+(l_before-l_after)
              end if
              !print the time of the recon and its location
              !also include the angle between the reconnecting filaments
              open(unit=61,file='./data/recon_location.log',position='append')
                write(61,*) t, 0.5*(f(i)%x+f(j)%x), acos(dot_val), (l_before-l_after)
              close(61)
              !now see if we setup reconnection distance array
              do k=1, n_recon_track
                if (full_recon_distance(k)%active.eqv..false.) then
                  full_recon_distance(k)%i=i !set i and 
                  full_recon_distance(k)%j=pari !j then say the array is 
                  full_recon_distance(k)%active=.true.!active and reset the
                  full_recon_distance(k)%counter=1 !counter
                  full_recon_distance(k)%angle=acos(dot_val) !set the angle
                  exit !we only set these particles in one of the slots
                end if
              end do 
            end if
            !set correct behind_infront
            f(parjb)%infront=pari
            f(pari)%behind=parjb
            f(i)%infront=j
            f(j)%behind=i
            if (recon_time_info) then
              !set reconnection times
              f(parjb)%t_recon(2)=f(parjb)%t_recon(1) !move last recon to slot 2
              f(parjb)%t_recon(1)=t !set current time
              f(pari)%t_recon(2)=f(pari)%t_recon(1) 
              f(pari)%t_recon(1)=t 
              f(i)%t_recon(2)=f(i)%t_recon(1) 
              f(i)%t_recon(1)=t 
              f(j)%t_recon(2)=f(j)%t_recon(1) 
              f(j)%t_recon(1)=t 
            end if
            !check the size of these new loops
            call loop_killer(pari) ; call loop_killer(i)
          end if
        end if 
      end if
    end do
  end subroutine
  !*********************************************************************
  !>reconnect filaments if they become too close very similar to default scheme
  !> but does not check line length before and after
  subroutine precon_non_dissipative
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
          !print*, 'reconnection'
          !print*, 'i',i, f(i)%infront, f(i)%behind
          !print*, 'j',j, f(j)%infront, f(j)%behind
          !reconnect the filaments
          recon_count=recon_count+1 !keep track of the total # of recons
          if (recon_info) then !more reconnection information
            call same_loop_test(i,j,same_loop) !general.mod
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
          !set reconnection times
          f(parjb)%t_recon(2)=f(parjb)%t_recon(1) !move last recon to slot 2
          f(parjb)%t_recon(1)=t !set current time
          f(pari)%t_recon(2)=f(pari)%t_recon(1) 
          f(pari)%t_recon(1)=t 
          f(i)%t_recon(2)=f(i)%t_recon(1) 
          f(i)%t_recon(1)=t 
          f(j)%t_recon(2)=f(j)%t_recon(1) 
          f(j)%t_recon(1)=t 
          !check the size of these new loops
          call loop_killer(pari) ; call loop_killer(i)
        end if 
      end if
    end do
  end subroutine
  !*********************************************************************
  !>reconnect filaments if they become too close using method of Kondaurova 
  !>et al. 2008 which requires that the filaments will pass through each other
  !>this is checked using a simple simultanious test
  subroutine precon_kondaurova
    use matrix
    implicit none
    real :: distr, min_distr !reconnection distances
    real :: dot_val, tangent1(3), tangent2(3) !used to determine if filaments parallel
    real :: l_before, l_after !to check line length before and after
    real :: kond_array(3,4), kond_soln(3) !for simultanious eqn solver
    integer :: kond_err !for simultanious eqn solver
    integer :: pari, parb, parii, parbb, parji, parjb !particles infront/behind
    integer :: par_recon !the particle we reconnect with
    integer :: i, j !we must do a double loop over all particles N^2
    logical :: kond_check, same_loop
    do i=1, pcount
      if (f(i)%infront==0) cycle !empty particle
      pari=f(i)%infront ; parb=f(i)%behind !find particle infront/behind
      parii=f(pari)%infront ; parbb=f(parb)%behind !find particle twice infront/behind
      !now we determine if we can reconnect
      if ((f(i)%closestd<2*delta).and.(f(i)%closestd>epsilon(1.))) then
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
        else
          !now we must check the line length before and after
          l_before=dist_gen(f(i)%x,f(i)%ghosti)+dist_gen(f(i)%x,f(i)%ghostb)+&
                   dist_gen(f(j)%x,f(j)%ghosti)+dist_gen(f(j)%x,f(j)%ghostb)
          l_after=dist_gen(f(i)%x,f(i)%ghostb)+dist_gen(f(i)%x,f(j)%x)+&
                   dist_gen(f(j)%x,f(j)%ghosti)+dist_gen(f(pari)%x,f(parjb)%x)
          if (l_after<=l_before) then
            !now we test wether the points will pass through each other
            !-------------column1--------------------
            kond_array(1,1)=f(i)%u(1)-f(j)%u(1) 
            kond_array(2,1)=f(i)%u(2)-f(j)%u(2)
            kond_array(3,1)=f(i)%u(3)-f(j)%u(3)
            !-------------column2--------------------
            kond_array(1,2)=f(i)%ghosti(1)-f(i)%x(1)
            kond_array(2,2)=f(i)%ghosti(2)-f(i)%x(2)
            kond_array(3,2)=f(i)%ghosti(3)-f(i)%x(3)
            !-------------column3--------------------
            kond_array(1,3)=-f(j)%ghosti(1)+f(j)%x(1)
            kond_array(2,3)=-f(j)%ghosti(2)+f(j)%x(2)
            kond_array(3,3)=-f(j)%ghosti(3)+f(j)%x(3)
            !-------------column4--------------------
            kond_array(1,4)=f(j)%x(1)-f(i)%x(1)
            kond_array(2,4)=f(j)%x(2)-f(i)%x(2)
            kond_array(3,4)=f(j)%x(3)-f(i)%x(3)
            call solve_linear_eqn(kond_array, kond_soln, 3, kond_err)           
            if ((kond_soln(1)>0.).and.(kond_soln(1)<dt).and. &
                (kond_soln(2)>0.).and.(kond_soln(2)<1.).and. &
                (kond_soln(3)>0.).and.(kond_soln(3)<1.)) then
              kond_check=.true.
            else
              kond_check=.false.
            end if
            if (kond_check) then
              !reconnect the filaments
              recon_count=recon_count+1 !keep track of the total # of recons
              if (recon_info) then !more reconnection information
                call same_loop_test(i,j,same_loop) !general.mod
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
              !set reconnection times
              f(parjb)%t_recon(2)=f(parjb)%t_recon(1) !move last recon to slot 2
              f(parjb)%t_recon(1)=t !set current time
              f(pari)%t_recon(2)=f(pari)%t_recon(1) 
              f(pari)%t_recon(1)=t 
              f(i)%t_recon(2)=f(i)%t_recon(1) 
              f(i)%t_recon(1)=t 
              f(j)%t_recon(2)=f(j)%t_recon(1) 
              f(j)%t_recon(1)=t 
              !check the size of these new loops
              call loop_killer(pari) ; call loop_killer(i)
            end if
          end if
        end if 
      end if
    end do
  end subroutine
  !**************************************************
  !>removes loops with less than 6 particles
  !>this is needed to ensure derivatives can be calculated correctly
  subroutine loop_killer(particle)
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
      if (mirror_bc) then
        !pinned particles will mess this up so exit routine
        !we have a separate test in mirror.mod
        if (f(next)%pinnedi.or.f(next)%pinnedb) return
      end if
    end do
    ! If loop is too small destroy
    if (counter<5) then
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
  !> small subroutine which allocates full_recon_distance 
  subroutine setup_recon_distance_array
    implicit none
    integer :: i
    write(*,*) 'length of reconnection distance array: ', size_recon_array
    write(*,*) 'number of pairs being tracked (simultaniously): ', n_recon_track
    allocate(full_recon_distance(n_recon_track))
    do i=1,n_recon_track
      allocate(full_recon_distance(i)%curvi(size_recon_array))
      allocate(full_recon_distance(i)%curvj(size_recon_array))
      allocate(full_recon_distance(i)%dist(size_recon_array))
    end do
    file_count_recon_array=1
  end subroutine
  !**************************************************
  !> called by run.f90 set the reconnection distance and print to file
  !> if the array is full 
  subroutine set_recon_dist
    implicit none
    character (len=40) :: print_file
    real :: dist
    integer :: i, k !for looping
    do k=1, n_recon_track !loop over all points we are tracking
      !skip k index if not active
      if (full_recon_distance(k)%active.eqv..false.) cycle
      !first check to see if either i or j has been removed
      if (f(full_recon_distance(k)%i)%infront==0) then
        !deactivate the array 
        full_recon_distance(k)%active=.false.
        !leave this comment in for now until completely satisfied
        print*, 'prematurely deacticvating', k
        return
      else if (f(full_recon_distance(k)%j)%infront==0) then
        !deactivate the array 
        full_recon_distance(k)%active=.false.
        !leave this comment in for now until completely satisfied
        print*, 'prematurely deacticvating', k
        return
      end if
      !now set the distances and curvatures
      dist=distf(full_recon_distance(k)%i,full_recon_distance(k)%j)
      if (dist>0.5*box_size) then
        !deactivate the array 
        full_recon_distance(k)%active=.false.
        !do an additional check if the distance array is 75% full we are better off
        !printing to file - note here we only write up to the current counter not the 
        !full array size
        write(unit=print_file,fmt="(a,i4.4,a)")"./data/recon_dist",file_count_recon_array,".log"
        open(unit=45,file=print_file,status='replace')
          write(45,*) t, '%time'
          write(45,*) dt, '%timestep'
          write(45,*) k, '% index of reconnection pair being tracked, (interest only)'
          write(45,*)  full_recon_distance(k)%counter, '%array size (reduced)'
          write(45,*) full_recon_distance(k)%angle, '%initial angle'
          do i=1, full_recon_distance(k)%counter
             write(45,*) full_recon_distance(k)%dist(i), full_recon_distance(k)%curvi(i), full_recon_distance(k)%curvj(i) 
          end do
        close(45)
        !increment file_count
        file_count_recon_array=file_count_recon_array+1
        return
      end if 
      full_recon_distance(k)%dist(full_recon_distance(k)%counter)=dist
      full_recon_distance(k)%curvi(full_recon_distance(k)%counter)=curvature(full_recon_distance(k)%i)
      full_recon_distance(k)%curvj(full_recon_distance(k)%counter)=curvature(full_recon_distance(k)%j)
      !increment counter
      full_recon_distance(k)%counter=full_recon_distance(k)%counter+1
      if (full_recon_distance(k)%counter==size_recon_array+1) then
        !deactivate the array 
        full_recon_distance(k)%active=.false.
        !print to file
        write(unit=print_file,fmt="(a,i4.4,a)")"./data/recon_dist",file_count_recon_array,".log"
        open(unit=45,file=print_file,status='replace')
          write(45,*) t, '%time'
          write(45,*) dt, '%timestep'
          write(45,*) k, '% index of reconnection pair being tracked, (interest only)'
          write(45,*)  size_recon_array, '%array size'
          write(45,*) full_recon_distance(k)%angle, '%initial angle'
          do i=1, size_recon_array
             write(45,*) full_recon_distance(k)%dist(i), full_recon_distance(k)%curvi(i), full_recon_distance(k)%curvj(i) 
          end do
        close(45)
        !increment file_count
        file_count_recon_array=file_count_recon_array+1
      end if
    end do
  end subroutine
end module
