function mesh_plot(filenumber)
close all
filename=sprintf('data/mesh%03d.dat',filenumber);
load data/dims.log;
msize=dims(3)
if (msize==0) 
  disp('mesh size is zero exiting script')
  return
end
fid=fopen(filename);
t=fread(fid,1,'float');
x=fread(fid,msize,'float');
unormx=fread(fid,msize^3,'float');
unormy=fread(fid,msize^3,'float');
unormz=fread(fid,msize^3,'float');
unorm_mrms=max(sqrt(unormx(:).^2+unormy(:).^2+unormz(:).^2));
ux=fread(fid,msize^3,'float');
uy=fread(fid,msize^3,'float');
uz=fread(fid,msize^3,'float');
unormx=reshape(unormx,msize,msize,msize);
unormy=reshape(unormy,msize,msize,msize);
unormz=reshape(unormz,msize,msize,msize);
ux=reshape(ux,msize,msize,msize);
uy=reshape(uy,msize,msize,msize);
uz=reshape(uz,msize,msize,msize);
%plot slices of field+isosurface
mesh_slices(x,ux,uy,uz,msize,'super')
if unorm_mrms>0.
  mesh_slices(x,unormx,unormy,unormz,msize,'normal')
end
%spectrum
mesh_spectrum(ux,uy,uz,msize,'super')
if unorm_mrms>0.
  mesh_spectrum(unormx,unormy,unormz,msize,'normal')
end
  
