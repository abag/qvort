function mesh_spec_fit(filenumber)
close all
filename=sprintf('data/mesh%03d.dat',filenumber);
load data/dims.log;
msize=dims(3)
if (msize==0) 
  disp('mesh size is zero exiting script')
  return
end
fid=fopen(filename);
t=fread(fid,1,'float64');
x=fread(fid,msize,'float64');
unormx=fread(fid,msize^3,'float64');
unormy=fread(fid,msize^3,'float64');
unormz=fread(fid,msize^3,'float64');
unorm_mrms=max(sqrt(unormx(:).^2+unormy(:).^2+unormz(:).^2));
ux=fread(fid,msize^3,'float64');
uy=fread(fid,msize^3,'float64');
uz=fread(fid,msize^3,'float64');
unormx=reshape(unormx,msize,msize,msize);
unormy=reshape(unormy,msize,msize,msize);
unormz=reshape(unormz,msize,msize,msize);
ux=reshape(ux,msize,msize,msize);
uy=reshape(uy,msize,msize,msize);
uz=reshape(uz,msize,msize,msize);
mesh_spectrum(ux,uy,uz,msize,'super',1)
if unorm_mrms>0.
  mesh_spectrum(unormx,unormy,unormz,msize,'normal',1)
end
  
