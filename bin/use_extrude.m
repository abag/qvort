function use_extrude(filenumber)
filename=sprintf('data/var%04d.log',filenumber);
%we set the dimensions of the box here
%this is overridden if we have periodic B.C.'s
box_size=.005 ;
%set options based on varargin
rough=0 ; linetrue=0 ; rainbow=0 ; dark=0 ; printit=0 ; overhead=0 ;
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
      break
  else
    varray(vcount,1)=x(next);
    varray(vcount,2)=y(next);
    varray(vcount,3)=z(next);
    next=round(f(next));
    vcount=vcount+1;
  end
end
%what shape do we extrude?
 q = linspace(0,2*pi,20);
 base = 0.5*[cos(q); sin(q)];   % Base curve is a circle (radius = 1)
 [X,Y,Z] = extrude(base,10000.*varray,1,[],1);
 figure; surf(X,Y,Z); axis equal;