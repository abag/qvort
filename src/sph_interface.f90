!>this module contains routines which interface between sph particles and the
!>filament
module sph_interface
  use cdata
  use general
  use sph
  implicit none
  contains
  !>use SPH particle velocity field
  subroutine SPH_f_interp()
    implicit none
    real :: u_interp(3) !interpolated velocity
    integer :: i, j !for looping
    real :: dist !distance between particles 
    do i=1, pcount
      if (f(i)%infront==0) cycle !empty particle
      if (f(i)%sph/=0) then 
        !particle is fixed to it's sph particle
        f(i)%x=s(f(i)%sph)%x
      else
        u_interp=0.
        do j=1, SPH_count
          dist=dist_gen(f(i)%x,s(j)%x) !distance between point and sph particle
          u_interp=u_interp+s(j)%u*s(j)%m*sph_W(dist,s(j)%h)/s(j)%rho
        end do
        !now timestep
        f(i)%x=f(i)%x+dt*u_interp
      end if
    end do
  end subroutine
  !**************************************************************************
  !>get divergence of SPH paricles at filament points
  subroutine SPH_f_divu(i,divu)
    implicit none
    real :: divu !divergence of u
    integer :: i,j !for looping
    real :: dist !distance between particles 
    if (f(i)%infront==0) return !empty particle
    if (f(i)%sph/=0) then
      !particle is fixed to it's sph particle
      divu=s(f(i)%sph)%divu
    else
      divu=0.
      do j=1, SPH_count
        dist=dist_gen(f(i)%x,s(j)%x) !distance between point and sph particle
        divu=divu+s(j)%divu*s(j)%m*sph_W(dist,s(j)%h)/s(j)%rho
      end do
    end if
  end subroutine
  !**************************************************************************
  !> is there an SPH particle our filament points can 'latch' onto, this is very
  !> similar to the reconnection alogorithm in line.mod. We loop over all particles
  !> and all vortex points and test to see if the distance between them is less than
  !> \f$\delta\f$, if so we associate the point with the sph particle.
  subroutine SPH_f_latch()
    implicit none
    integer :: i, j !for looping
    integer :: partner
    real :: dist !distance between particles 
    real :: min_dist
    do i=1, pcount
      if (f(i)%infront==0) cycle !empty particle
      if (f(i)%sph==0) then 
        !only do this for particles which are unassociated, i.e. f(i)%sph=0
        min_dist=10.
        do j=1, SPH_count
          dist=dist_gen(f(i)%x,s(j)%x) !distance between point and sph particle
          if (dist<min_dist) then
            min_dist=dist
            partner=j
          end if
        end do
        if (min_dist<delta) then
          f(i)%sph=partner !associate point i with sph particle j
        end if
      end if
    end do
  end subroutine
  !**************************************************************************
  !>interpolate particles between SPH particles using a gaussian process
  !>not needed at present but keep in the code as is a nice alternative to above
  subroutine SPH_f_GP()
    use matrix
    implicit none
    real :: u_interp(3) !interpolated velocity
    real :: disti, distb, t_herm !distances
    real :: mbehind, minfront
    real :: positions(4,3)
    real :: Kxstarx(4), Kxstarx_helper(4)
    real :: Kxx(4,4), Kxxinv(4,4) !gaussian process
    integer :: i, j, k !for looping
    integer :: minv_EF
    integer :: next
    integer :: infront, iinfront, behind, bbehind
    do i=1, pcount
      if (f(i)%infront==0) cycle !empty particle
      if (f(i)%sph/=0) then 
        !particle is fixed to it's sph particle
        f(i)%x=s(f(i)%sph)%x
      else
        !cubic hermite interpolation 
        !find particles infront, behind for interpolation
        !----------------infront-----------------
        next=f(i)%infront
        do j=1, pcount
          if (f(next)%sph==0) then
            next=f(next)%infront
          else
            infront=next
            exit
          end if
        end do
        !----------------iinfront-----------------
        next=f(infront)%infront
        do j=1, pcount
          if (f(next)%sph==0) then
            next=f(next)%infront
          else
            iinfront=next
            exit
          end if
        end do
        !----------------behind-----------------
        next=f(i)%behind
        do j=1, pcount
          if (f(next)%sph==0) then
            next=f(next)%behind
          else
            behind=next
            exit
          end if
        end do
        !----------------bbehind-----------------
        next=f(behind)%behind
        do j=1, pcount
          if (f(next)%sph==0) then
            next=f(next)%behind
          else
            bbehind=next
            exit
          end if
        end do
        positions(1,:)=f(bbehind)%x
        positions(2,:)=f(behind)%x
        positions(3,:)=f(infront)%x
        positions(4,:)=f(iinfront)%x
        !-----------------GP-----------------
        do j=1, 4
          do k=1,4
            Kxx(j,k)=covariance(dist_gen_sq(positions(j,:),positions(k,:)))
          end do
          Kxstarx(j)=covariance(dist_gen_sq(f(i)%x,positions(j,:)))
        end do
        !inverse matrix
        call findinv(Kxx, Kxxinv, 4, minv_EF) !matrix.mod
        if (minv_EF<0) then
          print*, '------------------matrix---------------'
          do j=1, 4
            print*, Kxx(j,:)
          end do
          print*, '------------------particles---------------'
          do j=1, 4
            print*, positions(j,:)
          end do
          print*, '------------------neighbours---------------'
          print*, 'i=',i,bbehind, behind, infront, iinfront
          call fatal_error('here','matrix printed above')
        end if 
        do j=1, 4
          Kxstarx_helper(j)=sum(Kxstarx*Kxxinv(j,:))
        end do
        f(i)%x=Kxstarx_helper(1)*positions(1,:)+&
               Kxstarx_helper(2)*positions(2,:)+&
               Kxstarx_helper(3)*positions(3,:)+&
               Kxstarx_helper(4)*positions(4,:)
      end if
    end do
  end subroutine
  !*************************************
  !>covariance function for Gaussian Process
  !>interpolation
  real function covariance(dist2)
    implicit none
    real, intent(IN) :: dist2
    covariance=exp(-dist2/((5*delta)**2))
  end function
end module
