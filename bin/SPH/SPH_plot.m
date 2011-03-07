function vortex_plot(filenumber)
filename=sprintf('data/SPH_par%04d.log',filenumber);
fid=fopen(filename);
%get the dimensions information from dims.log
dims=load('./data/dims.log');
if dims(4)==1
  %binary read
  t=fread(fid,1,'float64');
  number_of_particles=fread(fid,1,'int');
  x=fread(fid,number_of_particles,'float64');
  y=fread(fid,number_of_particles,'float64');
  z=fread(fid,number_of_particles,'float64');
  rho=fread(fid,number_of_particles,'float64');
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
    rho(j)=dummy_vect(4);
  end
end
range=ceil(10*max(rho)-10*min(rho))+1;
cmap=colormap(jet(range));
store_caxis=([min(rho) max(rho)]);
for i=1:number_of_particles
  plot3(x(i),y(i),z(i),'o','MarkerEdgeColor','k','MarkerFaceColor',cmap(1+ceil(10*rho(i)-10*min(rho)),:))
  hold on
end
axis([-dims(2)/2 dims(2)/2 -dims(2)/2 dims(2)/2 -dims(2)/2 dims(2)/2]); 
box on
caxis(store_caxis)
colorbar
set(gca,'FontSize',16)
