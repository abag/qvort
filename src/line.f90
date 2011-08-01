!>all routines which alter the geometry of the vortex filament/flux tube should
!>be contained in this module
module line
  use cdata
  use general
  use periodic
  use reconnection
  contains
  !>insert new points to maintain the resolution of the line,
  !>this will mean the point separation lies between \f$\delta\f$
  !>and \f$\delta/2\f$ to ensure the curvature is not affected by this the new point
  !>\f$ \mathbf{s}_{i'} \f$ is positioned between i, i+1 at position
  !>\f[ \mathbf{s}_{i'}=\frac{1}{2}(\mathbf{s}_i+\mathbf{s}_{i+1})+\left( \sqrt{R^2_{i'}
  !!-\frac{1}{4}\ell_{i+1}^2}-R_{i'} \right)\frac{\mathbf{s}_{i'}''}{|\mathbf{s}_{i'}''|},\f]
  !! where \f$R_{i'}=|\mathbf{s}_{i'}''|^{-1}\f$.
  subroutine pinsert
    implicit none    
    type(qvort), allocatable, dimension(:) :: tmp
    real :: disti, curv, f_ddot(3)
    integer :: par_new
    integer :: old_pcount, i
    !check that we have particles
    if (particles_only.eqv..false.) then
      if (count(mask=f(:)%infront>0)==0) then
        call fatal_error('line.mod:pinsert','vortex line length is 0, run over')
      end if
    else
     return !we don't need to run this routine if pcount=0
    end if
    old_pcount=pcount
    total_length=0. !zero this
    do i=1, old_pcount
      if (i>size(f)) then
        !bug check
        print*, 'I think there is a problem' ; exit
      end if
      if (f(i)%infront==0) cycle !empty particle
      if (mirror_bc.and.f(i)%pinnedi) cycle !pinned particle
      !get the distance between the particle and the one infront
      disti=dist_gen(f(i)%x,f(i)%ghosti) !general.f90
      total_length=total_length+disti !measure total length of filaments
      if (delta_adapt) then
        curv=curvature(i)
        f(i)%delta=1.+3.*exp(50.*(-delta*curv/2.))
      else
        f(i)%delta=1. !set the prefactor to 1
      end if
      if (disti>f(i)%delta*delta) then  
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
        curv=curv**(-1) !actually we want the inverse
        if (curv**2-0.25*disti**2>0.) then
          !could be negative, avoid this
          f(par_new)%x=0.5*(f(i)%x+f(i)%ghosti)+&
                       (sqrt(curv**2-0.25*disti**2)-curv)*curv*f_ddot
        else
          f(par_new)%x=0.5*(f(i)%x+f(i)%ghosti)
        end if
        !average the current velocity
        f(par_new)%u=0.5*(f(i)%u+f(f(i)%infront)%u) 
        !zero the older velocities
        f(par_new)%u1=0. ; f(par_new)%u2=0.
        !set correct infront and behinds & ghostzones
        f(par_new)%behind=i ; f(par_new)%infront=f(i)%infront
        call get_ghost_p(par_new,f(par_new)%ghosti, f(par_new)%ghostb,f(par_new)%ghostii, f(par_new)%ghostbb) !periodic.mod
        f(f(i)%infront)%behind=par_new
        call get_ghost_p(f(i)%infront,f(f(i)%infront)%ghosti, f(f(i)%infront)%ghostb, &
                                      f(f(i)%infront)%ghostii, f(f(i)%infront)%ghostbb) !periodic.mod         
        f(i)%infront=par_new
        call get_ghost_p(i,f(i)%ghosti, f(i)%ghostb,f(i)%ghostii, f(i)%ghostbb) !periodic.mod         
        !set local delta factor to be 1 incase we are adapting it
        f(par_new)%delta=1. 
      end if
    end do
    !calculate average separation of particles
    avg_sep=total_length/old_pcount
  end subroutine
  !*************************************************************************
  !>remove points along the filament if they are compressed to the point where
  !>the separation between the point i and i+2 is less than \f$\delta\f$
  !>if phonon emission is set to true in run.in then particles with high 
  !>curvature are removed, smoothing the loop to mimic dissipation at large k
  subroutine premove
    implicit none
    real :: distii
    integer :: infront, tinfront
    integer :: i
    logical :: do_remove 
    do i=1, pcount
      if (f(i)%infront==0) cycle !empty particle
      if (mirror_bc) then
        !do not test if you are pinned 
        if (f(i)%pinnedi) cycle 
        !or the particle infront is pinned
        if (f(f(i)%infront)%pinnedi) cycle 
      end if 
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
      if (distii<.499*delta*(f(i)%delta+f(f(i)%infront)%delta))then
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
      end if
    end do
  end subroutine
  !******************************************************************
  !>find the closest particle to i using N^2 operation this is
  !>done by looping over all particles and caling distf from general.mod
  !>\todo we need to do particles on the boundary as well (periodic)
  subroutine pclose
    implicit none
    integer :: i, j
    real :: dist
    do i=1, pcount
      if (f(i)%infront==0) cycle !empty particle
      f(i)%closestd=100. !arbitrarily high
      if (mirror_bc) then
        !do not test if you are pinned 
        if (f(i)%pinnedi.or.f(i)%pinnedb) cycle 
      end if 
      do j=1, pcount
        if (f(j)%infront==0) cycle !empty particle
        if ((i/=j).and.(f(i)%infront/=j).and.(f(i)%behind/=j)) then
        !the above line ensures we do not reconnect with particles infront/behind
          if (mirror_bc) then
            !do not test if j is pinned 
            if (f(j)%pinnedi.or.f(j)%pinnedb) cycle 
          end if 
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
  !>reconnect filaments if they become too close - this is a dummy routine which
  !> calls routines in the reconnection modul
  subroutine precon
    implicit none
    select case(recon_type)
      case('original')
        call precon_original !reconnection.mod
      case('dissipative')
        call precon_dissapitive !reconnection.mod
      case('non_dissipative')
        call precon_non_dissipative !reconnection.mod
      case('kondaurova')
        call precon_kondaurova !reconnection.mod
    end select
  end subroutine
end module
