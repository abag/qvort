module output
  !ALL THE OUTPUT/INPUT ROUTINES USED IN THE CODE
  use cdata
  contains
  !**********************************************************************
  subroutine print_dims()
    implicit none
    open(unit=77,file='./data/dims.log',status='replace')
      write(77,*) delta
      write(77,*) box_size
    close(77)  
  end subroutine
  !**********************************************************************
  subroutine print_info()
    !print information to screen
    implicit none
    if (itime==shots) then
      print*, 't---pcount--'
    end if
    write(*,'(f5.2,i5.4)') t, pcount
  end subroutine
  !**********************************************************************
  subroutine printf(filenumber)
    !print the f array as formatted data for use with gnuplot/matlab
    implicit none
    integer, intent(IN) :: filenumber
    character (len=40) :: print_file
    integer :: i
    write(unit=print_file,fmt="(a,i3.3,a)")"./data/var",filenumber,".log"
    open(unit=98,file=print_file,status='replace')
    write(98,*) t
    write(98,*) pcount
    do i=1, pcount
      write(98,*) f(i)%x(1:3), f(i)%infront
    end do
    close(98)
  end subroutine
  !**********************************************************************
  subroutine data_dump
    !store everything needed to restart the code
    implicit none
    open(unit=53,file="./data/var.dat",FORM='unformatted',status='replace')
      write(53) pcount
      write(53) itime
      write(53) t
      write(53) f
    close(53)
  end subroutine
  !**********************************************************************
  subroutine sdata_dump
    !store everything needed to restart the code - special edition :-)
    implicit none
    print*, 'dumping to special data file, current time is=', t
    open(unit=53,file="./data/special.dat",FORM='unformatted',status='replace')
      write(53) pcount
      write(53) itime
      write(53) t
      write(53) f
    close(53)
  end subroutine
  !**********************************************************************
  subroutine data_restore
    !restart the code - should this be in init.mod?
    implicit none
    integer :: dummy_itime 
    open(unit=63,file="./data/var.dat",FORM='unformatted')
      read(63) pcount
      read(63) itime
      read(63) t
      allocate(f(pcount))
      read(63) f
    close(63)
    nstart=itime+1
    print*, 'data read in from dump file at t=', t
  end subroutine
end module
