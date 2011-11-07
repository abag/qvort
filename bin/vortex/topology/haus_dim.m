function hausdorff=haus_dim(filenumber,option)
if nargin==1
  option='empty';
end
switch option
  case 'plot'
    disp('plotting loops with different colours')
  case 'empty'
  otherwise
    disp('incorrect option, aborting script and printing help:')
    return
end
filename=sprintf('data/var%04d.log',filenumber);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%get the dimensions information from dims.log
dims=load('./data/dims.log');
if dims(4)==1
  fid=fopen(filename);
  if fid<0
    disp('var file does not exist, exiting script')
    return
  end
  time=fread(fid,1,'float64');
  number_of_particles=fread(fid,1,'int');
  x=fread(fid,number_of_particles,'float64');
  y=fread(fid,number_of_particles,'float64');
  z=fread(fid,number_of_particles,'float64');
  f=fread(fid,number_of_particles,'int');
  u=fread(fid,number_of_particles,'float64');
  u2=fread(fid,number_of_particles,'float64');
else
  fid=fopen(filename);
  if fid<0
    disp('var file does not exist, exiting script')
    return
  end
  %read the time
  tline=fgetl(fid);
  dummy=textscan(tline, '%f');
  time=dummy{:};
  %how many particles
  tline=fgetl(fid);
  dummy=textscan(tline, '%d');
  number_of_particles=dummy{:};
  %get the particles%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  for j=1:number_of_particles
    tline=fgetl(fid);
    dummy=textscan(tline, '%f');
    dummy_vect=dummy{:};
    x(j)=dummy_vect(1);
    y(j)=dummy_vect(2);
    z(j)=dummy_vect(3);
    f(j)=dummy_vect(4);
    u(j)=dummy_vect(5);
    u2(j)=dummy_vect(6);
  end
  f=uint16(f);
end
zerocount=sum(f<=0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%find a starting point
for j=1:number_of_particles
  if f(j)~=0
    next=j;
    break
  end
end
%%%%%%%%%%%now do main loop%%%%%%%%%%%%%
line_count=0;
next_old=next;
counter(1:1000)=0;
haus_dim_count=0;
hausdorff=0;
for l=1:1000
  for i=1:number_of_particles
    line(l,i)=next;
    old_next=next;
    next=f(next);
    counter(l)=counter(l)+1;
    dummy_pos(counter(l),1)=x(old_next);
    dummy_pos(counter(l),2)=y(old_next);
    dummy_pos(counter(l),3)=z(old_next);
    if next==next_old
      break
    end
    if next==0
      break
    end
  end
  dummy_pos(counter(l)+1,1)=dummy_pos(1,1);  
  dummy_pos(counter(l)+1,2)=dummy_pos(1,2);
  dummy_pos(counter(l)+1,3)=dummy_pos(1,3);
  nbreakpoint=1;
  breakpoint(1)=0;
    for i=1:counter(l);
      %find distance between point and point infront
      dist=sqrt((dummy_pos(i+1,1)-dummy_pos(i,1))^2+(dummy_pos(i+1,2)-dummy_pos(i,2))^2+(dummy_pos(i+1,3)-dummy_pos(i,3))^2);
      if dist>0.5*dims(2)
        nbreakpoint=nbreakpoint+1;
        breakpoint(nbreakpoint)=i;
      end 
    end
    breakpoint(nbreakpoint+1)=counter(l);
  nbreakpoint;
  for i=1:nbreakpoint
 switch option
    case 'plot'
    plot3(dummy_pos(breakpoint(i)+1:breakpoint(i+1),1),dummy_pos(breakpoint(i)+1:breakpoint(i+1),2),dummy_pos(breakpoint(i)+1:breakpoint(i+1),3))
    hold on
end
    if (breakpoint(i+1)-breakpoint(i)-1)>10
    [n r]=boxcount(dummy_pos(breakpoint(i)+1:breakpoint(i+1),:));
    slope=polyfit(log(r),log(n),1);
    haus_dim_count=haus_dim_count+1;
    hausdorff=hausdorff-slope(1);
    end
  end
  clear breakpoint
  clear dummy_pos
  line_count=line_count+1;
  if sum(counter)<(number_of_particles-zerocount)
    for i=1:number_of_particles
      unique=true;
      for m=1:l
        for j=1:counter(m)
          if i==line(m,j)
            unique=false;
          end
          if f(i)==0
            unique=false;
          end
        end
      end
      if unique
        next=i;
        next_old=next;
      end
    end
  else
    break
  end
end
switch option
  case 'plot'
    hold off
    box on ; axis equal
    set(gca,'xtick',[]) ;  set(gca,'ytick',[]) ; set(gca,'ztick',[])
    rotate3d on
end
hausdorff=hausdorff/haus_dim_count;
