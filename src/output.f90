!>The main output routines from the code, whilst specifc modules may contain there own
!>output, the main routines for printing structures should be in here.
module output
  use cdata
  use tree
  contains
  !**********************************************************************
  !> print various dimensions/flags to file for matlab to read in for plotting
  subroutine print_dims()
    implicit none
    open(unit=77,file='./data/dims.log',status='replace')
      write(77,*) delta, '%resolution - \delta' 
      write(77,*) box_size, '%box size'
      write(77,*) mesh_size, '%mesh size'
      if (binary_print) then
        write(77,*) 1, '%binary printing'
      else
        write(77,*) 0, '%formatted printing'
      end if
      write(77,*) part_count, '%number of particles'
      write(77,*) quasi_pcount, '%number of quasi-particles'
      write(77,*) xdim_scaling_factor, '%length of x direction'
    close(77)
  end subroutine
  !**********************************************************************
  !>print information to screen/file
  subroutine print_info()
    implicit none
    call printf(itime/shots) !output.mod
    open(unit=78,file='data/ts.log',position='append')
    if (itime==shots) then
      write(*,'(a)') '--var--------t-------pcount-------recon----avg_d-----length&
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
    if (phonon_emission) then
      open(unit=78,file='data/phonon_count.log',position='append')
        write(78,*) phonon_count
      close(78)
   end if
    if (energy_inf) then
      open(unit=78,file='data/energy.log',position='append')
        write(78,*) energy
      close(78)
    end if
    if (recon_info) then
      open(unit=72,file='data/recon_extra_info.log',position='append')
        write(72,*) self_rcount, vv_rcount
      close(72)
    end if
    open(unit=79,file='data/curvature.log',position='append')
      write(79,*) kappa_bar, kappa_min, kappa_max
    close(79)
    if (topo_inf) then
      open(unit=78,file='./data/topology.log',position='append')
        write(78,*) t, linking_number, writhing_number
      close(78)
    end if
  end subroutine
  !**********************************************************************
  !>print the f (filament) array as (un)formatted data for use with gnuplot/matlab
  subroutine printf(filenumber)
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
        write(98) sqrt(f(:)%u_sup(1)**2+f(:)%u_sup(2)**2+f(:)%u_sup(3)**2)
      close(98)
    else
      write(unit=print_file,fmt="(a,i4.4,a)") "./data/var",filenumber,".log"
      open(unit=98,file=print_file,status='replace')
        write(98,*) t
        write(98,*) pcount
        do i=1, pcount
          write(98,*) f(i)%x(1:3), f(i)%infront, sqrt(f(i)%u(1)**2+f(i)%u(2)**2+f(i)%u(3)**2), &
                      sqrt(dot_product(f(i)%u_sup,f(i)%u_sup))
        end do
      close(98)
    end if
    if (simple_plots) then
      write(unit=print_file,fmt="(a,i4.4,a)")"./scripts/simple_plot.sh ",filenumber
      call system(print_file)
    end if
  end subroutine
  !**********************************************************************
  subroutine print_delta_adapt(filenumber)
    implicit none
    integer, intent(IN) :: filenumber
    character (len=40) :: print_file
    write(unit=print_file,fmt="(a,i4.4,a)")"./data/delta_adapt",filenumber,".dat"
    open(unit=98,file=print_file,status='replace',form='unformatted',access='stream')
      write(98) t
      write(98) pcount
      write(98) f(:)%delta
    close(98)
  end subroutine
  !**********************************************************************
  !>print the g (quasi particles) array as (un)formatted data for use with gnuplot/matlab
  subroutine printg(filenumber)
    implicit none
    integer, intent(IN) :: filenumber
    character (len=40) :: print_file
    integer :: i
    if (binary_print) then
      write(unit=print_file,fmt="(a,i4.4,a)")"./data/quasi_par",filenumber,".log"
      open(unit=98,file=print_file,status='replace',form='unformatted',access='stream')
        write(98) t
        write(98) quasi_pcount
        write(98) g(:)%x(1)
        write(98) g(:)%x(2)
        write(98) g(:)%x(3)
      close(98)
    else  
      write(unit=print_file,fmt="(a,i4.4,a)")"./data/quasi_par",filenumber,".log"
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
  !>print the p (particles) array as (un)formatted data for use with gnuplot/matlab
  subroutine printp(filenumber)
    implicit none
    integer, intent(IN) :: filenumber
    character (len=40) :: print_file
    integer :: i
    if (binary_print) then
      write(unit=print_file,fmt="(a,i4.4,a)")"./data/par",filenumber,".log"
      open(unit=98,file=print_file,status='replace',form='unformatted',access='stream')
        write(98) t
        write(98) part_count
        write(98) p(:)%x(1)
        write(98) p(:)%x(2)
        write(98) p(:)%x(3)
      close(98)
    else  
      write(unit=print_file,fmt="(a,i4.4,a)")"./data/par",filenumber,".log"
      open(unit=98,file=print_file,status='replace')
        write(98,*) t
        write(98,*) part_count
        do i=1, part_count
          write(98,*) p(i)%x(1:3)
        end do
      close(98)
    end if
  end subroutine
  !**********************************************************************
  !>store everything needed to restart the code
  !>\todo a few diagnostics are not dumped and resume from 0 please fix
  subroutine data_dump
    implicit none
    open(unit=53,file="./data/var.dat",FORM='unformatted',status='replace')
      write(53) pcount
      write(53) recon_count
      write(53) itime
      write(53) t
      write(53) f
      write(53) quasi_pcount
      if (quasi_pcount>0) write(53) g
      write(53) part_count
      if (part_count>0) write(53) p
    close(53)
  end subroutine
  !**********************************************************************
  !>store everything needed to restart the code - special edition :-)
  subroutine sdata_dump   
    implicit none
    write(*,*) 'dumping to special data file, current time is=', t
    open(unit=53,file="./data/special.dat",FORM='unformatted',status='replace')
      write(53) pcount
      write(53) recon_count
      write(53) itime
      write(53) t
      write(53) f
      write(53) quasi_pcount
      if (quasi_pcount>0) write(53) g
      write(53) part_count
      if (part_count>0) write(53) p
    close(53)
  end subroutine
  !**********************************************************************
  !>print the mesh to a binary file to be read in by matlab
  subroutine print_mesh(filenumber)
    implicit none
    integer, intent(IN) :: filenumber
    character (len=40) :: print_file
    real, allocatable :: vapor_array(:,:,:)
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
    !print the velocity field for vapor (vapor print set in run.in) 
    if (vapor_print) then
      allocate(vapor_array(mesh_size, mesh_size, mesh_size))
      vapor_array(:,:,:)=sqrt(mesh(:,:,:)%u_norm(1)**2+&
                              mesh(:,:,:)%u_norm(2)**2+&
                              mesh(:,:,:)%u_norm(3)**2)
      write(unit=print_file,fmt="(a,i3.3,a)")"./data/vap_norm",filenumber,".dat"
      open(unit=93,file=print_file,form='unformatted',status='replace',access='stream')
        write(93) vapor_array
      close(93)
      vapor_array(:,:,:)=sqrt(mesh(:,:,:)%u_sup(1)**2+&
                              mesh(:,:,:)%u_sup(2)**2+&
                              mesh(:,:,:)%u_sup(3)**2)
      write(unit=print_file,fmt="(a,i3.3,a)")"./data/vap_sup",filenumber,".dat"
      open(unit=93,file=print_file,form='unformatted',status='replace',access='stream')
        write(93) vapor_array
      close(93)
      deallocate(vapor_array) 
    end if
  end subroutine
  !**********************************************************************
  !>print the velocity  array as unformatted data for use matlab
  subroutine print_velocity(filenumber)
    implicit none
    integer, intent(IN) :: filenumber
    character (len=40) :: print_file
    integer :: i
    write(unit=print_file,fmt="(a,i4.4,a)")"./data/uu",filenumber,".dat"
    open(unit=98,file=print_file,status='replace',form='unformatted',access='stream')
      write(98) t
      write(98) pcount
      write(98) f(:)%u_sup(1)
      write(98) f(:)%u_sup(2)
      write(98) f(:)%u_sup(3)
    close(98)
    if (vel_print_extra) then
      write(unit=print_file,fmt="(a,i4.4,a)")"./data/uu_full",filenumber,".dat"
      open(unit=98,file=print_file,status='replace',form='unformatted',access='stream')
        write(98) t
        write(98) pcount
        write(98) f(:)%u(1)
        write(98) f(:)%u(2)
        write(98) f(:)%u(3)
      close(98)
      write(unit=print_file,fmt="(a,i4.4,a)")"./data/uu_mf",filenumber,".dat"
      open(unit=98,file=print_file,status='replace',form='unformatted',access='stream')
        write(98) t
        write(98) pcount
        write(98) f(:)%u_mf(1)
        write(98) f(:)%u_mf(2)
        write(98) f(:)%u_mf(3)
      close(98)
    end if
  end subroutine
  !**************************************************
  subroutine banner_print()
    implicit none
    write(*,*) "                                 _   " 
    write(*,*) "            __ ___   _____  _ __| |_ "
    write(*,*) "           / _` \ \ / / _ \| '__| __|"
    write(*,*) "          | (_| |\ V / (_) | |  | |_ "
    write(*,*) "           \__, | \_/ \___/|_|   \__|"
    write(*,*) "              |_|                    "
    write(*,*) "                                     " 
  end subroutine
end module
