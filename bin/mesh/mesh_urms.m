function [t urms nurms]=mesh_urms(filenumber)
filename=sprintf('./data/mesh%03d.dat',filenumber);
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
t=fread(fid,1,'float64');
x=fread(fid,msize,'float64');
unormx=fread(fid,msize^3,'float64');
unormy=fread(fid,msize^3,'float64');
unormz=fread(fid,msize^3,'float64');
nurms=sqrt(sum(unormx(:).^2+unormy(:).^2+unormz(:).^2)/msize^3);
ux=fread(fid,msize^3,'float64');
uy=fread(fid,msize^3,'float64');
uz=fread(fid,msize^3,'float64');
urms=sqrt(sum(ux(:).^2+uy(:).^2+uz(:).^2)/msize^3);
  
