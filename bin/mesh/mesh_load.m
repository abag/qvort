function mesh_load(filenumber)
global ux uy uz msize unormx unormy unormz unorm_mrms u_mrms x t
load data/dims.log;
msize=dims(3);
if (msize==0) 
  disp('mesh size is zero exiting script')
  return
end
filename=sprintf('./data/mesh%03d.dat',filenumber);
fid=fopen(filename);
if fid<0
  disp('mesh file does not exist, exiting script')
  return
end
disp(sprintf('mesh size is: %04d',msize))
t=fread(fid,1,'float64');
x=fread(fid,msize,'float64');
unormx=fread(fid,msize^3,'float64');
unormy=fread(fid,msize^3,'float64');
unormz=fread(fid,msize^3,'float64');
unorm_mrms=max(sqrt(unormx(:).^2+unormy(:).^2+unormz(:).^2));
ux=fread(fid,msize^3,'float64');
uy=fread(fid,msize^3,'float64');
uz=fread(fid,msize^3,'float64');
u_mrms=max(sqrt(ux(:).^2+uy(:).^2+uz(:).^2));
unormx=reshape(unormx,msize,msize,msize);
unormy=reshape(unormy,msize,msize,msize);
unormz=reshape(unormz,msize,msize,msize);
ux=reshape(ux,msize,msize,msize);
uy=reshape(uy,msize,msize,msize);
uz=reshape(uz,msize,msize,msize);