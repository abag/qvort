module output
  !MAIN OUTPUT ROUTINES USED IN THE CODE
  use cdata
  use tree
  contains
  !**********************************************************************
  subroutine print_dims()
    implicit none
    open(unit=77,file='./data/dims.log',status='replace')
      write(77,*) delta
      write(77,*) box_size
      write(77,*) mesh_size
      if (binary_print) then
        write(77,*) 1
      else
        write(77,*) 0
      end if
    close(77)
  end subroutine
  !**********************************************************************
  subroutine print_info()
    !print information to screen/file
    implicit none
    open(unit=78,file='data/ts.log',position='append')
    if (itime==shots) then
      write(*,*) '--var--------t-------pcount-------recon----avg_d-----length&
                  --------maxu---------maxdu-----num eval----curv------removed'
      write(78,*) '%--var--------t-------pcount-------recon----avg_d-----length&
                   --------maxu---------maxdu-----num eval----curv------removed'
    end if
    write(*,'(i6.4,f13.7,i10.1,i13.1,f7.4,f13.6,f13.5,f13.5,f10.2,f10.2,i13.1)') &
itime/shots,t,count(mask=f(:)%infront>0),recon_count,avg_sep/delta,&
total_length,maxu,maxdu,real(eval_counter)/count(mask=f(:)%infront>0),kappa_bar,&
remove_count
    write(78,'(i6.4,f13.7,i10.1,i13.1,f7.4,f13.6,f13.5,f13.5,f10.2,f10.2,i13.1)') &
itime/shots,t,count(mask=f(:)%infront>0),recon_count,avg_sep/delta,&
total_length,maxu,maxdu,real(eval_counter)/count(mask=f(:)%infront>0),kappa_bar,&
remove_count
    close(78)
    open(unit=78,file='data/energy.log',position='append')
      write(78,*) energy
    close(78)
    open(unit=79,file='data/curvature.log',position='append')
      write(79,*) kappa_bar, kappa_min, kappa_max
    close(79)
  end subroutine
  !**********************************************************************
  subroutine printf(filenumber)
    !print the f array as (un)formatted data for use with gnuplot/matlab
    implicit none
    integer, intent(IN) :: filenumber
    character (len=40) :: print_file
    integer :: i
    if (filenumber==10000) call warning_message('output.mod','run out of filenumbers to print var to')
    if (binary_print) then
      write(unit=print_file,fmt="(a,i4.4,a)")"./data/var",filenumber,".log"
      open(unit=98,file=print_file,status='replace',form='unformatted',access='stream')
        write(98) t
        write(98) pcount
        write(98) f(:)%x(1)
        write(98) f(:)%x(2)
        write(98) f(:)%x(3)
        write(98) f(:)%infront
        write(98) sqrt(f(:)%u(1)**2+f(:)%u(2)**2+f(:)%u(3)**2)
      close(98)
    else
      write(unit=print_file,fmt="(a,i4.4,a)")"./data/var",filenumber,".log"
      open(unit=98,file=print_file,status='replace')
        write(98,*) t
        write(98,*) pcount
        do i=1, pcount
          write(98,*) f(i)%x(1:3), f(i)%infront, sqrt(f(i)%u(1)**2+f(i)%u(2)**2+f(i)%u(3)**2)
        end do
      close(98)
    end if
  end subroutine
  !**********************************************************************
  subroutine printg(filenumber)
    !print the g (particles) array as (un)formatted data for use with gnuplot/matlab
    implicit none
    integer, intent(IN) :: filenumber
    character (len=40) :: print_file
    integer :: i
    if (binary_print) then
      write(unit=print_file,fmt="(a,i4.4,a)")"./data/par",filenumber,".log"
      open(unit=98,file=print_file,status='replace',form='unformatted',access='stream')
        write(98) t
        write(98) quasi_pcount
        write(98) g(:)%x(1)
        write(98) g(:)%x(2)
        write(98) g(:)%x(3)
      close(98)
    else  
      write(unit=print_file,fmt="(a,i4.4,a)")"./data/par",filenumber,".log"
      open(unit=98,file=print_file,status='replace')
        write(98,*) t
        write(98,*) quasi_pcount
        do i=1, quasi_pcount
          write(98,*) g(i)%x(1:3)
        end do
      close(98)
    end if
  end subroutine
  !**********************************************************************
  subroutine data_dump
    !store everything needed to restart the code
    implicit none
    open(unit=53,file="./data/var.dat",FORM='unformatted',status='replace')
      write(53) pcount
      write(53) recon_count
      write(53) itime
      write(53) t
      write(53) f
    close(53)
  end subroutine
  !**********************************************************************
  subroutine sdata_dump
    !store everything needed to restart the code - special edition :-)
    implicit none
    write(*,*) 'dumping to special data file, current time is=', t
    open(unit=53,file="./data/special.dat",FORM='unformatted',status='replace')
      write(53) pcount
      write(53) itime
      write(53) t
      write(53) f
    close(53)
  end subroutine
  !**********************************************************************
  subroutine print_mesh(filenumber)
    !print the mesh to a binary file
    implicit none
    integer, intent(IN) :: filenumber
    character (len=40) :: print_file
    real, allocatable :: vapor_array(:,:,:)
    logical :: vapor=.true.
    if (mesh_size==0) return
    write(unit=print_file,fmt="(a,i3.3,a)")"./data/mesh",filenumber,".dat"
    open(unit=92,file=print_file,form='unformatted',status='replace',access='stream')
      write(92) t
      write(92) mesh(mesh_size/2,mesh_size/2,1:mesh_size)%x(1)
      write(92) mesh(1:mesh_size,1:mesh_size,1:mesh_size)%u_norm(1)
      write(92) mesh(1:mesh_size,1:mesh_size,1:mesh_size)%u_norm(2)
      write(92) mesh(1:mesh_size,1:mesh_size,1:mesh_size)%u_norm(3)
      write(92) mesh(1:mesh_size,1:mesh_size,1:mesh_size)%u_sup(1)
      write(92) mesh(1:mesh_size,1:mesh_size,1:mesh_size)%u_sup(2)
      write(92) mesh(1:mesh_size,1:mesh_size,1:mesh_size)%u_sup(3)
    close(92)
    !print just the velocity field for vapor 
    if (vapor) then
      allocate(vapor_array(mesh_size, mesh_size, mesh_size))
      vapor_array(:,:,:)=sqrt(mesh(:,:,:)%u_sup(1)**2+&
                              mesh(:,:,:)%u_sup(2)**2+&
                              mesh(:,:,:)%u_sup(3)**2)
      write(unit=print_file,fmt="(a,i3.3,a)")"./data/vap_mesh",filenumber,".dat"
      open(unit=93,file=print_file,form='unformatted',status='replace',access='stream')
        write(93) vapor_array
      close(93)
      deallocate(vapor_array) 
    end if
  end subroutine
  !**********************************************************************
  subroutine print_velocity(filenumber)
    !print the velocity  array as unformatted data for use matlab
    implicit none
    integer, intent(IN) :: filenumber
    character (len=40) :: print_file
    integer :: i
    if (vel_print) then !set in run.in false by default
      write(unit=print_file,fmt="(a,i4.4,a)")"./data/uu",filenumber,".dat"
      open(unit=98,file=print_file,status='replace',form='unformatted',access='stream')
        write(98) t
        write(98) pcount
        write(98) f(:)%u(1)
        write(98) f(:)%u(2)
        write(98) f(:)%u(3)
      close(98)
    end if
  end subroutine
end module
