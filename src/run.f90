program run
  use cdata
  use initial
  use output
  use periodic
  use timestep
  use line
  implicit none
  logical :: restart
  call init_random_seed !cdata.mod
  !read in parameters
  call read_run_file !cdata.mod
  !initial conditions-can we restart?
  inquire(file="./data/var.dat", exist=restart)
  if (restart) then
    call data_restore !output.mod
  else
    call init_setup !initial.mod
  end if
  !print dimensions info for matlab plotting
  call print_dims !output.mod
  !begin time loop
  do itime=nstart, nsteps
    !---------------------main operations--------------------------
    call ghostp !periodic.mod
    call pmotion !timestep.mod
    call pinsert !line.mod
    if (mod(itime, recon_shots)==0) then
      call premove !line.mod
    end if
    if(periodic_bc) call enforce_periodic !periodic.mod
    !--------------now do all data output--------------------------
    if (mod(itime, shots)==0) then
      !print to file
      call printf(itime/shots) !output.mod
      !store data to a binary file for reload
      call data_dump !output.mod
      call print_info !output.mod
    end if
    !special data dump
    if ((itime/=0).and.(itime==int_special_dump)) call sdata_dump
    !--------------------------------------------------------------
  end do
!still to do:-
!biot savart
!analysis
end program
