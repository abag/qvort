function smoothed_field(filenumber,pview)
close all
if (nargin==1)
  pview=0
end
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
if (pview==1)
  disp('printing to raw mesh for paraview/vapor')
  %savevtkvector(wx,wy,wz,'pout_smooth_vector.vtk');
  %savevtk(sqrt(wx.^2+wy.^2+wz.^2),'pout_smooth.vtk');
  writevtk(sqrt(wx.^2+wy.^2+wz.^2),'pout_smooth.vtk');
  return
end
%return
%plot slices of field+isosurface
mesh_slices(x,wx,wy,wz,msize,'smoothed')
%current
[jx jy jz]=curl(wx,wy,wz);
mesh_slices(x,jx,jy,jz,msize,'current')

