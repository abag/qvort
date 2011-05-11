function mesh_plot(filenumber,varargin)
close all
optargin = size(varargin,2);
%set options based on varargin
slice=0 ; iso=0 ; spect=0 ; struct=0 ; print=0 ; 
for i=1:optargin
  switch cell2str(varargin(i))
    case 'slice'
      slice=1;
    case 'iso'
      iso=1;
    case 'spect'
      spect=1;
    case 'struct'
      struct=1;
    case 'print'
      print=1;
      disp('printing to file')
  end
end
filename=sprintf('data/mesh%03d.dat',filenumber);
load data/dims.log;
msize=dims(3);
if (msize==0) 
  disp('mesh size is zero exiting script')
  return
end
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
%plot slices of field
if slice==1
  if u_mrms>0.
    mesh_slices(x,ux,uy,uz,msize,'super')
  end
  if unorm_mrms>0.
    mesh_slices(x,unormx,unormy,unormz,msize,'normal')
  end
end
%plot isosurface
if iso==1
  if u_mrms>0.
    mesh_iso(x,ux,uy,uz,msize,'super')
  end
  if unorm_mrms>0.
    mesh_iso(x,unormx,unormy,unormz,msize,'normal')
  end
end
%spectrum
if u_mrms>0.
  if spect==1
    mesh_spectrum(ux,uy,uz,msize,'super',0)
  end
  if struct==1
    %mesh_structure_func(x,ux,uy,uz,msize,'super')
  end
end
if unorm_mrms>0.
  if spect==1
    mesh_spectrum(unormx,unormy,unormz,msize,'normal',0)
  end
  if struct==1
    mesh_structure_func(x,unormx,unormy,unormz,msize,'normal')
  end
end
  
