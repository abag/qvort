!>all reconnection routines are contained in this module this then feeds into line.mod
module reconnection
  use cdata
  use general
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
  !>reconnect filaments if they become too close - default scheme
  subroutine precon_original
    implicit none
    real :: distr, min_distr !reconnection distances
    real :: dot_val, tangent1(3), tangent2(3) !used to determine if filaments parallel
    real :: l_before, l_after !to check line length before and after
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
  !*********************************************************************
  !>reconnect filaments if they become too close using method of Kondaurova 
  !>et al. 2008 which requires that the filaments will pass through each other
  !>this is checked using a simple simultanious test
  subroutine precon_kondaurova
    implicit none
    real :: distr, min_distr !reconnection distances
    real :: dot_val, tangent1(3), tangent2(3) !used to determine if filaments parallel
    real :: l_before, l_after !to check line length before and after
    real :: kond_a, kond_b, kond_c, kond_d !for simultanious eqn solver
    real :: kond_e, kond_f, kond_g, kond_h !for simultanious eqn solver
    real :: kond_x, kond_y !for simultanious eqn solver
    real :: kond_check !test simultanious solution for x/y is valid in z
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
        else
          !now we must check the line length before and after
          l_before=dist_gen(f(i)%x,f(i)%ghosti)+dist_gen(f(i)%x,f(i)%ghostb)+&
                   dist_gen(f(j)%x,f(j)%ghosti)+dist_gen(f(j)%x,f(j)%ghostb)
          l_after=dist_gen(f(i)%x,f(i)%ghostb)+dist_gen(f(i)%x,f(j)%x)+&
                   dist_gen(f(j)%x,f(j)%ghosti)+dist_gen(f(pari)%x,f(parjb)%x)
          if (l_after<=l_before) then
            !now we test wether the points will pass through each other
            kond_a=f(i)%x(1)+dt*f(i)%u(1)
            kond_b=f(i)%ghosti(1)-f(i)%x(1)
            kond_c=f(j)%x(1)+dt*f(j)%u(1)
            kond_d=f(j)%ghosti(1)-f(j)%x(1)
            kond_e=f(i)%x(2)+dt*f(i)%u(2)
            kond_f=f(i)%ghosti(2)-f(i)%x(2)
            kond_g=f(j)%x(2)+dt*f(j)%u(2)
            kond_h=f(j)%ghosti(2)-f(j)%x(2)
            !above was setting variables now we solve simultanious eqn
            kond_y=(kond_e+(kond_f/kond_b)*(kond_c-kond_a)-kond_g)/(kond_h-(kond_f*kond_d/kond_b))
            kond_x=(kond_c-kond_a+kond_d*kond_y)/kond_b
            !now we check these in the final equation for z
            kond_check=f(i)%x(3)+dt*f(i)%u(3)+kond_x*(f(i)%ghosti(3)-f(i)%x(3))-&
                      (f(j)%x(3)+dt*f(j)%u(3)+kond_y*(f(j)%ghosti(3)-f(j)%x(3)))
            if (isnan(kond_check)) return !can happen if points lie on a plane
            if (kond_check<1E-4) then !THIS TOLERANCE NEEDS TO BE CHECKED
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
end module
