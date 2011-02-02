%function compressible
close all
filename='data/compressible_mesh.dat';
msize=128;
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
ux=fread(fid,msize^3,'float64');
uy=fread(fid,msize^3,'float64');
uz=fread(fid,msize^3,'float64');
fclose(fid)
ux=reshape(ux,msize,msize,msize);
uy=reshape(uy,msize,msize,msize);
uz=reshape(uz,msize,msize,msize);
unorm=max(sqrt(ux.^2+uy.^2+uz.^2));
%plot slices of field+isosurface
mesh_slices(x,ux,uy,uz,msize,'compressible')

  
