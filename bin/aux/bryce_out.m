function use_extrude(filenumber)
filename=sprintf('data/var%04d.log',filenumber);
%we set the dimensions of the box here
%this is overridden if we have periodic B.C.'s
box_size=.005 ;
%set options based on varargin
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
  end
  f=uint16(f);
end
%now create vectors to plot%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vcount=1;
next=1; % start at the first particle
for j=1:number_of_particles
  if round(f(j))==0
      f(j)
      break
  else
    varray(vcount,1)=x(next);
    varray(vcount,2)=y(next);
    varray(vcount,3)=z(next);
    varray(vcount,4)=.00015;
    next=round(f(next));
    vcount=vcount+1;
  end
end
varray(vcount,1:3)=varray(1,1:3);
varray(vcount,4)=varray(1,4);
vcount1=1:length(varray)
vcount2=1:0.05:length(varray);
varray2(:,1)=spline(vcount1,varray(:,1),vcount2);
varray2(:,2)=spline(vcount1,varray(:,2),vcount2);
varray2(:,3)=spline(vcount1,varray(:,3),vcount2);
varray2(:,4)=spline(vcount1,varray(:,4),vcount2);
saveobjtube('vortex_bryce.obj',varray2(:,1),varray2(:,2),varray2(:,3),varray2(:,4))
%saveobjtube('vortex_bryce.obj',varray2(:,1),varray2(:,2),varray2(:,3))