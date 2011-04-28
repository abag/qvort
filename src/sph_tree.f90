!>all the tree routines used by the SPH code are contained within this module,
module sph_tree
   use cdata
   use general
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
   type SPH_node
      real :: posx,posy,posz 
      real :: width 
      integer :: pcount 
      type (SPH_node), pointer :: fbl, fbr !front (y) bottom left, right (child)
      type (SPH_node), pointer :: ftl, ftr !front (y) top left, right (child)
      type (SPH_node), pointer :: bbl, bbr !back (y) bottom left, right (child)
      type (SPH_node), pointer :: btl, btr !back (y) top left, right (child)
      type (smooth_particle), allocatable :: parray(:) 
      real :: centx, centy, centz 
      real :: mass !total circulation vector
   end type SPH_node
   !>the SPH particle tree
   type (SPH_node), pointer :: stree 
  contains
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
  subroutine SPH_new_tree(stree,posx,posy,posz,width, loc_pcount,parray)
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
    stree%centz=0. ; stree%mass=0.
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
        !also account for the circulation
        stree%mass=stree%mass+parray(i)%m
      end if
    end do
    !take the average of the positions
    if (counter/=0) then
     stree%centx=stree%centx/counter
     stree%centy=stree%centy/counter
     stree%centz=stree%centz/counter
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
end module
