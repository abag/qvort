!> topological aspects of vortex filament
!> included in diagnostic.mod
module topology
  use cdata
  use general
  implicit none
  contains
    !>get the linking number of the current filament
    subroutine get_linking_number
      implicit none
      integer :: i, j
      integer :: crossings, crossings_check
      real, dimension(3) :: v1, v2, u1, u2
      integer :: crossing_detail
      logical :: same_loop
      crossings=0
      crossings_check=0
      do i=1, pcount
        if (f(i)%infront==0) cycle
        do j=1, pcount
          if (i==j) cycle
          if (f(j)%infront==0) cycle
          call same_loop_test(i,j,same_loop)
          if (same_loop) cycle
          ! we now define the four points that we
          ! test for a crossing projection in xy plane
          call xy_projection(v1,v2,u1,u2,i,j)
          !now we call our crossing algorithm
          call crossing_routine(v1,v2,u1,u2,crossing_detail)
          crossings=crossings+crossing_detail
          call zx_projection(v1,v2,u1,u2,i,j)
          !now we call our crossing algorithm
          call crossing_routine(v1,v2,u1,u2,crossing_detail)
          crossings_check=crossings_check+crossing_detail
        end do
      end do
      linking_number=crossings/2
    end subroutine
    !*************************************************************************
    !>get the writhing number of the current filament, to do this the an
    !!integral must be taken over the solid angle, we do this here in 40
    !!iterations which seems enough for convergence.
    subroutine get_writhing_number
      implicit none
      real, dimension(3) :: v1, v2, u1, u2
      real :: rot_angle
      integer :: i, j, k
      integer :: infronti, infrontj
      integer :: crossings
      integer :: crossing_detail
      logical :: same_loop
      crossings=0
      do k=1, 20 !20 projections is sufficient
        do i=1, pcount
          if (f(i)%infront==0) cycle
          do j=1, pcount
            if (i==j) cycle
            if (f(j)%infront==0) cycle
            call same_loop_test(i,j,same_loop)
            if (same_loop.eqv..false.) cycle
            infronti=f(i)%infront ; infrontj=f(j)%infront
            ! we now define the four points that we
            ! test for a crossing, however for writhing
            ! number we must integrate over solid angle
            rot_angle=(k-1)*2*pi/20.
            v1(1)=f(i)%x(1)
            v1(2)=f(i)%x(2)*cos(rot_angle)+f(i)%x(3)*sin(rot_angle)
            v1(3)=f(i)%x(2)*(-sin(rot_angle))+f(i)%x(3)*cos(rot_angle)
            v2(1)=f(infronti)%x(1)
            v2(2)=f(infronti)%x(2)*cos(rot_angle)+f(infronti)%x(3)*sin(rot_angle)
            v2(3)=f(infronti)%x(2)*(-sin(rot_angle))+f(infronti)%x(3)*cos(rot_angle)
            u1(1)=f(j)%x(1)
            u1(2)=f(j)%x(2)*cos(rot_angle)+f(j)%x(3)*sin(rot_angle)
            u1(3)=f(j)%x(2)*(-sin(rot_angle))+f(j)%x(3)*cos(rot_angle)
            u2(1)=f(infrontj)%x(1)
            u2(2)=f(infrontj)%x(2)*cos(rot_angle)+f(infrontj)%x(3)*sin(rot_angle)
            u2(3)=f(infrontj)%x(2)*(-sin(rot_angle))+f(infrontj)%x(3)*cos(rot_angle)
            !now we call our crossing algorithm
            call crossing_routine(v1,v2,u1,u2,crossing_detail)
            crossings=crossings+crossing_detail
          end do
        end do
      end do
      do k=1, 20
        do i=1, pcount
          if (f(i)%infront==0) cycle
          do j=1, pcount
            if (i==j) cycle
            if (f(j)%infront==0) cycle
            call same_loop_test(i,j,same_loop)
            if (same_loop.eqv..false.) cycle
            infronti=f(i)%infront ; infrontj=f(j)%infront
            ! we now define the four points that we
            ! test for a crossing, however for writhing
            ! number we must integrate over solid angle
            rot_angle=(k-1)*2*pi/20.
            v1(1)=f(i)%x(1)*cos(rot_angle)+f(i)%x(3)*(-sin(rot_angle))
            v1(2)=f(i)%x(2)
            v1(3)=f(i)%x(1)*sin(rot_angle)+f(i)%x(3)*cos(rot_angle)
            v2(1)=f(infronti)%x(1)*cos(rot_angle)+f(infronti)%x(3)*(-sin(rot_angle))
            v2(2)=f(infronti)%x(2)
            v2(3)=f(infronti)%x(1)*sin(rot_angle)+f(infronti)%x(3)*cos(rot_angle)
            u1(1)=f(j)%x(1)*cos(rot_angle)+f(j)%x(3)*(-sin(rot_angle))
            u1(2)=f(j)%x(2)
            u1(3)=f(j)%x(1)*sin(rot_angle)+f(j)%x(3)*cos(rot_angle)
            u2(1)=f(infrontj)%x(1)*cos(rot_angle)+f(infrontj)%x(3)*(-sin(rot_angle))
            u2(2)=f(infrontj)%x(2)
            u2(3)=f(infrontj)%x(1)*sin(rot_angle)+f(infrontj)%x(3)*cos(rot_angle)
            !now we call our crossing algorithm
            call crossing_routine(v1,v2,u1,u2,crossing_detail)
            crossings=crossings+crossing_detail
          end do
        end do
      end do
      writhing_number=real(crossings)/40
    end subroutine
    !*************************************************************************
    !>test if a set of 2D vectors cross and if so what is the direction of the 
    !>crossing.
    subroutine crossing_routine(v1,v2,u1,u2,crossing_detail)
      implicit none
      real, dimension(3) :: v1, v2, u1, u2
      real, dimension(2) :: v1_rot,v2_rot,u1_rot,u2_rot, intersect
      real :: a1, a2, b1, b2
      real :: lambda, mu
      real :: z1, z2
      real :: crossp
      integer :: crossing_detail
      crossing_detail=0
      if ((v1(1)==v2(1)).or.(u1(1)==u2(1))) then
        !print*, 'vertical lines', i, j
        u1_rot(1)=u1(1)*cos(0.7)-u1(2)*sin(0.7)
        u1_rot(2)=u1(1)*sin(0.7)+u1(2)*cos(0.7)
        
        u2_rot(1)=u2(1)*cos(0.7)-u2(2)*sin(0.7)
        u2_rot(2)=u2(1)*sin(0.7)+u2(2)*cos(0.7)
        
        v1_rot(1)=v1(1)*cos(0.7)-v1(2)*sin(0.7)
        v1_rot(2)=v1(1)*sin(0.7)+v1(2)*cos(0.7)
        
        v2_rot(1)=v2(1)*cos(0.7)-v2(2)*sin(0.7)
        v2_rot(2)=v2(1)*sin(0.7)+v2(2)*cos(0.7)
        
        u1(1)=u1_rot(1)
        u1(2)=u1_rot(2)
        u2(1)=u2_rot(1) 
        u2(2)=u2_rot(2)
        v1(1)=v1_rot(1)
        v1(2)=v1_rot(2)
        v2(1)=v2_rot(1)
        v2(2)=v2_rot(2)
      end if
      !construct a's and b's
      b1=(v2(2)-v1(2))/(v2(1)-v1(1))
      b2=(u2(2)-u1(2))/(u2(1)-u1(1))
      a1=v1(2)-b1*v1(1)
      a2=u1(2)-b2*u1(1)

      if (b1==b2) then
        !print*, 'parallel'
      else
        intersect(1)=-(a1-a2)/(b1-b2)
        intersect(2)=a1+b1*intersect(1)
        if (((v1(1)-intersect(1))*(intersect(1)-v2(1)))>=0.and. &
           ((u1(1)-intersect(1))*(intersect(1)-u2(1)))>=0.and. &
           ((v1(2)-intersect(2))*(intersect(2)-v2(2)))>=0.and. &
           ((u1(2)-intersect(2))*(intersect(2)-u2(2)))>=0) then
          ! determine wether an underpass or overpass
          lambda=(intersect(1)-v1(1))/(v2(1)-v1(1))
          mu=(intersect(1)-u1(1))/(u2(1)-u1(1))
          z1=v1(3)+lambda*(v2(3)-v1(3))
          z2=u1(3)+mu*(u2(3)-u1(3))
          if (z1<z2) then
            !underpass
            !need to take cross product to determine wether to left or to right
            crossp=(v2(1)-v1(1))*(u2(2)-u1(2)) &
                         -(v2(2)-v1(2))*(u2(1)-u1(1))
            if (crossp<0)  crossing_detail=1
            if (crossp>0)  crossing_detail=-1
          end if
        end if         
      end if         
    end subroutine
    !*************************************************************************
    !>take 4 points (in 3D space) - in pairs and apply a slight rotation
    !> then take xy projection
    subroutine xy_projection(v1,v2,u1,u2,i,j)
      implicit none
      real, dimension(3) :: v1, v2, u1, u2
      real :: rot_angle=0.8!rotate the data slightly for numerical accuracy
      integer, intent(IN) :: i, j
      integer:: infronti, infrontj
      infronti=f(i)%infront ; infrontj=f(j)%infront
      v1(1)=f(i)%x(1)
      v1(2)=f(i)%x(2)*cos(rot_angle)+f(i)%x(3)*sin(rot_angle)
      v1(3)=f(i)%x(2)*(-sin(rot_angle))+f(i)%x(3)*cos(rot_angle)
      v2(1)=f(infronti)%x(1)
      v2(2)=f(infronti)%x(2)*cos(rot_angle)+f(infronti)%x(3)*sin(rot_angle)
      v2(3)=f(infronti)%x(2)*(-sin(rot_angle))+f(infronti)%x(3)*cos(rot_angle)
 
      u1(1)=f(j)%x(1)
      u1(2)=f(j)%x(2)*cos(rot_angle)+f(j)%x(3)*sin(rot_angle)
      u1(3)=f(j)%x(2)*(-sin(rot_angle))+f(j)%x(3)*cos(rot_angle)
      u2(1)=f(infrontj)%x(1)
      u2(2)=f(infrontj)%x(2)*cos(rot_angle)+f(infrontj)%x(3)*sin(rot_angle)
      u2(3)=f(infrontj)%x(2)*(-sin(rot_angle))+f(infrontj)%x(3)*cos(rot_angle)
    end subroutine
    !*************************************************************************
    !>take 4 points (in 3D space) and project onto zx plane
    subroutine zx_projection(v1,v2,u1,u2,i,j)
      implicit none
      real, dimension(3) :: v1, v2, u1, u2
      integer, intent(IN) :: i, j
      integer:: infronti, infrontj
      infronti=f(i)%infront ; infrontj=f(j)%infront
      v1(3)=f(i)%x(1)
      v1(1)=f(i)%x(2)
      v1(2)=f(i)%x(3)
      v2(3)=f(infronti)%x(1)
      v2(1)=f(infronti)%x(2)
      v2(2)=f(infronti)%x(3)
      u1(3)=f(j)%x(1)
      u1(1)=f(j)%x(2)
      u1(2)=f(j)%x(3)
      u2(3)=f(infrontj)%x(1)
      u2(1)=f(infrontj)%x(2)
      u2(2)=f(infrontj)%x(3)
    end subroutine
    !*************************************************************************
    !>take 4 points (in 3D space) and project onto yz plane
    subroutine yz_projection(v1,v2,u1,u2,i,j)
      implicit none
      real, dimension(3) :: v1, v2, u1, u2
      integer, intent(IN) :: i, j
      integer:: infronti, infrontj
      infronti=f(i)%infront ; infrontj=f(j)%infront
      v1(2)=f(i)%x(1)
      v1(3)=f(i)%x(2)
      v1(1)=f(i)%x(3)
      v2(2)=f(infronti)%x(1)
      v2(3)=f(infronti)%x(2)
      v2(1)=f(infronti)%x(3)
      u1(2)=f(j)%x(1)
      u1(3)=f(j)%x(2)
      u1(1)=f(j)%x(3)
      u2(2)=f(infrontj)%x(1)
      u2(3)=f(infrontj)%x(2)
      u2(1)=f(infrontj)%x(3)
    end subroutine
endmodule topology
