function norm_init
close all
filename='data/norm_init_mesh.dat';
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
ux=fread(fid,msize^3,'float64');
uy=fread(fid,msize^3,'float64');
uz=fread(fid,msize^3,'float64');
fclose(fid)
ux=reshape(ux,msize,msize,msize);
uy=reshape(uy,msize,msize,msize);
uz=reshape(uz,msize,msize,msize);
unorm=max(sqrt(ux.^2+uy.^2+uz.^2));
figure('Name','u(x)')
  plot(squeeze(ux(32,32,:)))
%plot slices of field+isosurface
mesh_slices(x,ux,uy,uz,msize,'initial normal fluid')
%plot slices of divergence of field
figure('Name','div(u)')
div = divergence(ux,uy,uz);
pcolor(squeeze(div(32,:,:))) ; shading interp
colorbar
%plot slices of curl of field
figure('Name','curl(u)')
curlz = curl(ux,uy,uz);
pcolor(squeeze(curlz(32,:,:))) ; shading interp
figure
isosurface(sqrt(ux.^2+uy.^2+uz.^2))
