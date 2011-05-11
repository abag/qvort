!>all the tree routines used by the SPH code are contained within this module,
module sph_tree
   use cdata
   use general
   use sph_kernel
   !>the tree structure, uses pointers
   !>@param posx @param posy @param posz bottem left corner of box
   !>@param width the width of the mesh 
   !>@param fbl child cell front bottom left, @param fbr child cell front bottom right
   !>@param ftl child cell front top left @param ftr child cell front top right
   !>@param bbl child cell back bottom left, @param bbr child cell back bottom right
   !>@param btl child cell back top left @param btr child cell back top right
   !>@param all particles in the cell
   !>@param centx @param centy @param centz the center of 'mass' of the cell
   !>@param mass the total mass of the cell
   !>@param h the average smoothing length of the cell
   type SPH_node
      integer :: i !sph particle only set if pcount is 1
      real :: posx,posy,posz 
      real :: width 
      integer :: pcount 
      type (SPH_node), pointer :: fbl, fbr !front (y) bottom left, right (child)
      type (SPH_node), pointer :: ftl, ftr !front (y) top left, right (child)
      type (SPH_node), pointer :: bbl, bbr !back (y) bottom left, right (child)
      type (SPH_node), pointer :: btl, btr !back (y) top left, right (child)
      type (smooth_particle), allocatable :: parray(:) 
      real :: centx, centy, centz 
      real :: mass !total mass 
      real :: h
   end type SPH_node
   !>the SPH particle tree
   type (SPH_node), pointer :: stree
   integer, private :: ncounter
  contains
  !**************************************************************************
  !>An array to set the main node of the tree and then call the recurisive 
  !>subroutine build_tree this is called in the main program
  subroutine create_SPH_tree()
    allocate(stree)
    allocate(stree%parray(SPH_count))
    !set up the main node
    stree%posx=-box_size/2. ; stree%posy=-box_size/2. ; stree%posz=-box_size/2.
    stree%width=box_size
    stree%pcount=SPH_count
    stree%parray=s
    stree%centx=sum(s%x(1))/SPH_count
    stree%centy=sum(s%x(2))/SPH_count
    stree%centz=sum(s%x(3))/SPH_count
    stree%mass=sum(s%m)
    stree%h=sum(s%h)/SPH_count
    !now nullify its children - for safety
    nullify(stree%fbl, stree%fbr, stree%ftl, stree%ftr)
    nullify(stree%bbl, stree%bbr, stree%btl, stree%btr)
    call SPH_build_tree(stree) !the recursive routine which builds the tree
  end subroutine
  !**************************************************************************
  !>called by above, build the tree, continue to subdivide until a cell is empty
  !>or contains only one particle - octree structure used each new cell is created
  !>by the routine below new_tree
  recursive subroutine SPH_build_tree(stree)
    implicit none
    type(SPH_node), pointer :: stree
    character (len=40) :: print_file
    if (stree%pcount>1) then
      !divide up into 8 children
      !front, bottom, left (x,y,z)
      call SPH_new_tree(stree%fbl, stree%posx, stree%posy, stree%posz, &
                    stree%width/2., stree%pcount, stree%parray)
      !front, bottom, right (x+w,y,z)
      call SPH_new_tree(stree%fbr, stree%posx+(stree%width/2), stree%posy, stree%posz, &
                    stree%width/2., stree%pcount, stree%parray)
      !front, top, left (x,y,z+w)
      call SPH_new_tree(stree%ftl, stree%posx, stree%posy,stree%posz+(stree%width/2), &
                    stree%width/2., stree%pcount, stree%parray)
      !front, top, right (x+w,y,z+w)
      call SPH_new_tree(stree%ftr, stree%posx+(stree%width/2), stree%posy,stree%posz+(stree%width/2), &
                    stree%width/2., stree%pcount, stree%parray)
      !back, bottom, left (x,y+w,z)
      call SPH_new_tree(stree%bbl, stree%posx, stree%posy+(stree%width/2), stree%posz, &
                    stree%width/2., stree%pcount, stree%parray)
      !back, bottom, right (x+w,y+w,z)
      call SPH_new_tree(stree%bbr, stree%posx+(stree%width/2), stree%posy+(stree%width/2), stree%posz, &
                    stree%width/2., stree%pcount, stree%parray)
      !back, top, left (x,y+w,z+w)
      call SPH_new_tree(stree%btl, stree%posx, stree%posy+(stree%width/2), stree%posz+(stree%width/2), &
                    stree%width/2., stree%pcount, stree%parray)
      !back, top, right (x+w,y+w,z+w)
      call SPH_new_tree(stree%btr, stree%posx+(stree%width/2), stree%posy+(stree%width/2), stree%posz+(stree%width/2), &
                    stree%width/2., stree%pcount, stree%parray)
    else
      !print to file to test the algorithm works
      if ((mod(itime,shots)==0).and.tree_print) then
        write(unit=print_file,fmt="(a,i3.3,a)")"./data/tree_mesh",itime/shots,".log"
        open(unit=97,file=print_file,action='write', position='append')
          write(97,*) stree%posx, stree%posx+stree%width, stree%posy, &
                     stree%posy+stree%width, stree%posz, stree%posz+stree%width, stree%pcount
        close(97)
      end if
    end if
  end subroutine
  !********************************************************************************
  !>allocate a new cell, populate with particles and calculate centre of mass, 
  !>total circulation vector and magnetic field vector
  subroutine SPH_new_tree(stree,posx,posy,posz,width,loc_pcount,parray)
    !create a new tree
    implicit none
    integer, intent(IN) :: loc_pcount
    real, intent(IN) :: posx, posy, posz, width
    type(smooth_particle) :: parray(loc_pcount)
    type(SPH_node), pointer :: stree
    integer :: i, counter
    allocate(stree)
    allocate(stree%parray(loc_pcount))
    !set up this node of the tree with the input info
    stree%posx=posx ; stree%posy=posy ; stree%posz=posz
    stree%width=width ; stree%centx=0. ; stree%centy=0.
    stree%centz=0. ; stree%mass=0. ; stree%h=0.
    counter=0 !0 the counter
    do i=1, loc_pcount
      if((parray(i)%x(1)>=posx).and.(parray(i)%x(1)<(posx+width)).and. &
         (parray(i)%x(2)>=posy).and.(parray(i)%x(2)<(posy+width)).and. &
         (parray(i)%x(3)>=posz).and.(parray(i)%x(3)<(posz+width))) then
        counter=counter+1
        stree%parray(counter)=parray(i)
        !add the position of the particle to the position of the mesh
        stree%centx=stree%centx+parray(i)%x(1)
        stree%centy=stree%centy+parray(i)%x(2)
        stree%centz=stree%centz+parray(i)%x(3)
        !also account for the mass
        stree%mass=stree%mass+parray(i)%m
        !and the smoothing length
        stree%h=stree%h+parray(i)%h
      end if
    end do
    !take the average of the positions
    if (counter/=0) then
     stree%centx=stree%centx/counter
     stree%centy=stree%centy/counter
     stree%centz=stree%centz/counter
     !and average the smoothing length
     stree%h=stree%h/counter
    end if
    stree%pcount=counter
    !nullify the children
    nullify(stree%fbl, stree%fbr, stree%ftl, stree%ftr)
    nullify(stree%bbl, stree%bbr, stree%btl, stree%btr)
    call SPH_build_tree(stree) !recall the routine that builds the tree
  end subroutine
  !********************************************************************************
  !>empty out the tree structure to avoid a memory leak - repeatedly calls dead_tree
  recursive subroutine SPH_empty_tree(stree)
    implicit none
    type (SPH_node), pointer :: stree
    if(stree%pcount>1) then
      !kill all it's 8 children
      call SPH_dead_tree(stree%fbl)
      call SPH_dead_tree(stree%fbr)
      call SPH_dead_tree(stree%ftl)
      call SPH_dead_tree(stree%ftr)
      call SPH_dead_tree(stree%bbl)
      call SPH_dead_tree(stree%bbr)
      call SPH_dead_tree(stree%btl)
      call SPH_dead_tree(stree%btr)
    end if
  end subroutine
  !********************************************************************************
  !>deallocate a tree
  subroutine SPH_dead_tree(stree)
    !chop down that tree (metaphorically) 
    implicit none
    type (SPH_node), pointer :: stree
    call SPH_empty_tree(stree) !kill its children (if it has any)
    deallocate(stree%parray)
    deallocate(stree)
    nullify(stree)
  end subroutine
  !********************************************************************************
  !>find the closest particle using the tree structure, theoretically should be NlogN
  !>repeatedly calls SPH_tree_neighbour
  subroutine SPH_neighbour_list
    !loop over all the particles and find the closest one to i using the tree mesh
    implicit none
    integer :: i
    do i=1, SPH_count
      ncounter=0
      nullify(s_NN(i)%list)
      allocate(s_NN(i)%list)
      nullify(s_NN(i)%list%next)
      s_NN(i)%current => s_NN(i)%list
      call SPH_tree_neighbour(i,stree)
      s(i)%ncount=ncounter
    end do
  end subroutine
  !********************************************************************************
  !>recurisve routine used by SPH_neighbour list to find the neighbours of particle to i
  recursive subroutine SPH_tree_neighbour(i,stree)
    implicit none
    integer, intent(IN) :: i
    type (SPH_node), pointer :: stree
    real :: dist !distance
     if (stree%pcount==0) return !empty box no use
     !work out distances opening angles etc. here
     dist=sqrt((s(i)%x(1)-stree%centx)**2+(s(i)%x(2)-stree%centy)**2+(s(i)%x(3)-stree%centz)**2)
     if (dist-sqrt(3.)*stree%width<2*s(i)%h) then
       !if the box is empty exit the routine
       if (stree%pcount==0) then 
         return
       elseif (stree%pcount==1) then
         !if box only contains one particle then we add it to the neighbour list
         ncounter=ncounter+1
         if (ncounter==1) then
           s_NN(i)%current%i=stree%parray(1)%i
         else
           allocate(s_NN(i)%current%next)
           nullify(s_NN(i)%current%next%next)
           s_NN(i)%current%next%i=stree%parray(1)%i
           s_NN(i)%current => s_NN(i)%current%next
         end if
       else
         !open the box up and use the child cells
         call SPH_tree_neighbour(i,stree%fbl) ; call SPH_tree_neighbour(i,stree%fbr)
         call SPH_tree_neighbour(i,stree%ftl) ; call SPH_tree_neighbour(i,stree%ftr)
         call SPH_tree_neighbour(i,stree%bbl) ; call SPH_tree_neighbour(i,stree%bbr)
         call SPH_tree_neighbour(i,stree%btl) ; call SPH_tree_neighbour(i,stree%btr)
       end if
     end if
   end subroutine
   !************************************************************************
   !>a more general form of above, where a position is an input
   recursive subroutine SPH_gravity_tree_walk(i,stree,a)
     implicit none
     integer, intent(IN) :: i !the SPH particle
     type (SPH_node), pointer :: stree !the tree
     real :: a(3) !the acceleration due to G at particle i
     real :: vect(3) !helper vector
     real :: dist, theta !helper variables
     integer :: j=0 !the particle in the tree mesh

     if (stree%pcount==0) return !empty box no use
     !work out distances opening angles etc. here
     dist=sqrt((s(i)%x(1)-stree%centx)**2+&
               (s(i)%x(2)-stree%centy)**2+&
               (s(i)%x(3)-stree%centz)**2)
     theta=stree%width/dist !the angle to test below
 
     if (stree%pcount==1.or.theta<SPH_theta) then
       !check the distance is not 0
       if (dist<epsilon(0.)) return 
       !use the contribution of this cell
       vect(1)=s(i)%x(1)-stree%centx
       vect(2)=s(i)%x(2)-stree%centy
       vect(3)=s(i)%x(3)-stree%centz
       
       a=a-vect*SPH_G*stree%mass*sph_phi(dist,s(i)%h)/(2.*dist)
       a=a-vect*SPH_G*stree%mass*sph_phi(dist,stree%h)/(2.*dist)
     else
       !open the box up and use the child cells
       call SPH_gravity_tree_walk(i,stree%fbl,a) ; call SPH_gravity_tree_walk(i,stree%fbr,a)
       call SPH_gravity_tree_walk(i,stree%ftl,a) ; call SPH_gravity_tree_walk(i,stree%ftr,a)
       call SPH_gravity_tree_walk(i,stree%bbl,a) ; call SPH_gravity_tree_walk(i,stree%bbr,a)
       call SPH_gravity_tree_walk(i,stree%btl,a) ; call SPH_gravity_tree_walk(i,stree%btr,a)
     end if 
   end subroutine
end module
