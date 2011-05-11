function mesh_plot(filenumber)
close all
filename=sprintf('data/SPH_mesh%03d.dat',filenumber);
load data/dims.log;
msize=dims(8)
if (msize==0) 
  disp('mesh size is zero exiting script')
  return
end
fid=fopen(filename);
if fid<0
  disp('mesh file does not exist, exiting script')
  return
end
t=fread(fid,1,'float64');
x=fread(fid,msize,'float64');
rho=fread(fid,msize^3,'float64');
rho=reshape(rho,msize,msize,msize);
pcolor(squeeze(rho(1:msize,1:msize,msize/2)))
colorbar
%plot slices of field+isosurface
  
