function vortex_plot(filenumber)
filename=sprintf('data/par%03d.log',filenumber);
fid=fopen(filename);
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
%get the dimensions information from dims.log
dims=load('./data/dims.log');
plot3(x,y,z,'o','MarkerEdgeColor','k','MarkerFaceColor','r')
if (dims(2)>0.)
  axis([-dims(2)/2 dims(2)/2 -dims(2)/2 dims(2)/2 -dims(2)/2 dims(2)/2]); 
  box on
else
  axis([-0.1 0.1 -0.1 0.1 -0.1 0.1]);      
end
