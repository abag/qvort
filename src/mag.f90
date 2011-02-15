module mag
  !MAGNETIC FIELD 
  use cdata
  use general
  use smoothing
  use tree
  contains
 !**********************************************************************
  subroutine mag_tension(i,u)
    implicit none
    integer,intent(IN) :: i !particle
    real, intent(OUT) :: u(3) !tension velocity
    real, dimension(3) :: B,Bpx,Bmx,Bpy,Bmy,Bpz,Bmz,J
    logical, parameter :: fullJxB=.false.
    if (fullJxB) then
      B=0. ; call tree_smooth(f(i)%x,vtree,(/0.,0.,0./),B)
      Bpx=0. ; call tree_smooth(f(i)%x+(/delta,0.,0./),vtree,(/0.,0.,0./),Bpx)
      Bmx=0. ; call tree_smooth(f(i)%x-(/delta,0.,0./),vtree,(/0.,0.,0./),Bmx)
      Bpy=0. ; call tree_smooth(f(i)%x+(/0.,delta,0./),vtree,(/0.,0.,0./),Bpy)
      Bmy=0. ; call tree_smooth(f(i)%x-(/0.,delta,0./),vtree,(/0.,0.,0./),Bmy)
      Bpz=0. ; call tree_smooth(f(i)%x+(/0.,0.,delta/),vtree,(/0.,0.,0./),Bpz)
      Bmz=0. ; call tree_smooth(f(i)%x-(/0.,0.,delta/),vtree,(/0.,0.,0./),Bmz)
      J(1)=Bpy(3)-Bmy(3)-Bpz(2)+Bmz(2)
      J(2)=Bpz(1)-Bmz(1)-Bpx(3)+Bmx(3)
      J(3)=Bpx(2)-Bmx(2)-Bpy(1)+Bmy(1)
      J=J/(2*delta)
      u=cross_product(J,B)
    else
      call get_deriv_2(i,u)
      u=u*(f(i)%B**2)
    end if
  end subroutine
  !**********************************************************************
  subroutine B_smooth()
    implicit none
    real, allocatable :: B(:)
    integer :: i
    allocate(B(pcount))
    do i=1, pcount
      if (f(i)%infront==0) then
        B(i)=0. !empty particle
      else
        !locally smooth the field
        B(i)=0.25*(f(f(i)%behind)%B+2*f(i)%B+f(f(i)%infront)%B)
      end if
    end do
    f(:)%B=B
  end subroutine
  !**********************************************************************
  subroutine B_ts
    implicit none
    open(unit=71,file='data/B_ts.log',position='append')
      if (itime==shots) then
        write(71,*) '%--------t--------Bmax---------Bmin------------Brms-------'
      end if
      write(71,'(f13.7,f13.7,f13.7,f13.7)') t, maxval(f(:)%B), minval(f(:)%B), Brms
    close(71)
  end subroutine
  !**********************************************************************
  subroutine Bstretch(filenumber)
    implicit none
    integer, intent(IN) :: filenumber
    character (len=40) :: print_file
    integer :: i
    write(unit=print_file,fmt="(a,i4.4,a)")"./data/Bstretch",filenumber,".log"
    open(unit=98,file=print_file,status='replace')
    do i=1,pcount
      if (f(i)%infront==0) cycle
      write(98,*) f(i)%B
    end do
    close(98)
  end subroutine
end module
