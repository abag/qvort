function mesh_plot(filenumber,varargin)
close all
optargin = size(varargin,2);
%set options based on varargin
slice=0 ; current_slice=0. ; iso=0 ; para=0 ; render=0 ; 
for i=1:optargin
  switch cell2str(varargin(i))
    case 'slice'
      slice=1;
    case 'current_slice'
      current_slice=1;
    case 'iso'
      iso=1;
    case 'render'
      render=1;
    case 'para'
      para=1;
  end
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
if para==1
  disp('printing to raw mesh for paraview/vapor')
  %savevtkvector(wx,wy,wz,'pout_smooth_vector.vtk');
  %savevtk(sqrt(wx.^2+wy.^2+wz.^2),'pout_smooth.vtk');
  writevtk(sqrt(wx.^2+wy.^2+wz.^2),'pout_smooth.vtk');
end
%plot slices of field+isosurface
if slice==1
  mesh_slices(x,wx,wy,wz,msize,'smoothed')
end
if current_slice==1
  [jx jy jz]=curl(wx,wy,wz);  
  mesh_slices(x,jx,jy,jz,msize,'current')
end
if iso==1
  mesh_iso(x,wx,wy,wz,msize,'smoothed')
end
if render==1
  VolumeRender(sqrt(wx.^2+wy.^2+wz.^2)); 
end
rotate3d on

