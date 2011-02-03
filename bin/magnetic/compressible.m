function compressible
close all
filename='data/compressible_mesh.dat';
load data/nfm_dims.log;
msize=nfm_dims(1)
if (msize==0) 
  disp('mesh size is zero exiting script')
  return
end
fid=fopen(filename);
if fid<0
  disp('mesh file does not exist, exiting script')
  return
end
x=fread(fid,msize,'float64');
v=fread(fid,msize^3,'float64');
div=fread(fid,msize^3,'float64');
ux=fread(fid,msize^3,'float64');
uy=fread(fid,msize^3,'float64');
uz=fread(fid,msize^3,'float64');
fclose(fid)
v=reshape(v,msize,msize,msize);
div=reshape(div,msize,msize,msize);
ux=reshape(ux,msize,msize,msize);
uy=reshape(uy,msize,msize,msize);
uz=reshape(uz,msize,msize,msize);
unorm=max(sqrt(ux.^2+uy.^2+uz.^2));
%plot slices of field+isosurface
mesh_slices(x,ux,uy,uz,msize,'forcing: vel field')
scalar_slices(x,v,msize,'forcing: scalar field')
scalar_slices(x,div,msize,'forcing: divergence')

  
