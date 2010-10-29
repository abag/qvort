function vortex_plot(filenumber)
filename=sprintf('data/par%04d.log',filenumber);
fid=fopen(filename);
%get the dimensions information from dims.log
dims=load('./data/dims.log');
if dims(4)==1
  t=fread(fid,1,'float64');
  number_of_particles=fread(fid,1,'int');
  x=fread(fid,number_of_particles,'float64');
  y=fread(fid,number_of_particles,'float64');
  z=fread(fid,number_of_particles,'float64');
else 
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
  end
end
plot3(x,y,z,'o','MarkerEdgeColor','k','MarkerFaceColor','r')
if (dims(2)>0.)
  axis([-dims(2)/2 dims(2)/2 -dims(2)/2 dims(2)/2 -dims(2)/2 dims(2)/2]); 
  box on
else
  axis([-0.1 0.1 -0.1 0.1 -0.1 0.1]);      
end
