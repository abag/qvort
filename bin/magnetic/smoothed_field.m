function smoothed_field(filenumber)
close all
filename=sprintf('data/smoothed_field%03d.dat',filenumber);
load data/sm_dims.log;
msize=sm_dims(1)
fid=fopen(filename);
if fid<0
  disp('file does not exist, exiting script')
  return
end
t=fread(fid,1,'float64');
x=fread(fid,msize,'float64');
wx=fread(fid,msize^3,'float64');
wy=fread(fid,msize^3,'float64');
wz=fread(fid,msize^3,'float64');
wx=reshape(wx,msize,msize,msize);
wy=reshape(wy,msize,msize,msize);
wz=reshape(wz,msize,msize,msize);
%plot slices of field+isosurface
mesh_slices(x,wx,wy,wz,msize,'smoothed')
  
