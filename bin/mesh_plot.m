function mesh_plot(filenumber)
close all
filename=sprintf('data/mesh%03d.dat',filenumber);
load data/dims.log;
msize=dims(3);
fid=fopen(filename);
t=fread(fid,1,'float');
x=fread(fid,msize,'float');
ux=fread(fid,msize^3,'float');
uy=fread(fid,msize^3,'float');
uz=fread(fid,msize^3,'float');
ux=reshape(ux,msize,msize,msize);
uy=reshape(uy,msize,msize,msize);
uz=reshape(uz,msize,msize,msize);
%plot slices of field+isosurface
mesh_slices(x,ux,uy,uz,msize)
%spectrum
mesh_spectrum(ux,uy,uz,msize)
  
