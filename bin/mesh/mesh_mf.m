function mesh_mf(filenumber)
load data/dims.log;
msize=dims(3);
if (msize==0) 
  disp('mesh size is zero exiting script')
  return
end
filename=sprintf('./data/mesh_mf%03d.dat',filenumber);
fid=fopen(filename);
if fid<0
  disp('mesh file does not exist, exiting script')
  return
end
disp(sprintf('mesh size is: %04d',msize))
t=fread(fid,1,'float64');
x=fread(fid,msize,'float64');
mf1x=fread(fid,msize^3,'float64');
mf1y=fread(fid,msize^3,'float64');
mf1z=fread(fid,msize^3,'float64');
mf2x=fread(fid,msize^3,'float64');
mf2y=fread(fid,msize^3,'float64');
mf2z=fread(fid,msize^3,'float64');
mf1x=reshape(mf1x,msize,msize,msize);
mf1y=reshape(mf1y,msize,msize,msize);
mf1z=reshape(mf1z,msize,msize,msize);
mf2x=reshape(mf2x,msize,msize,msize);
mf2y=reshape(mf2y,msize,msize,msize);
mf2z=reshape(mf2z,msize,msize,msize);
max(max(max(mf1y)))
mf1x=2.79E-5*mf1x+(1.96E-5-9.97E-4*2.36E-2)*mf2x;
mf1y=2.79E-5*mf1y+(1.96E-5-9.97E-4*2.36E-2)*mf2y;
mf1z=2.79E-5*mf1z+(1.96E-5-9.97E-4*2.36E-2)*mf2z;
max(sqrt(mf1x.^2+mf1y.^2+mf1z.^2))
index = find(sqrt(mf1x.^2+mf1y.^2+mf1z.^2) > 10);
%mf1x(index)=0.;
%mf1y(index)=0.;
%mf1z(index)=0.;
%ksdensity(reshape(sqrt(mf1x.^2+mf1y.^2+mf1z.^2),msize^3,1))
%return
%figure
mesh_iso(x,mf1x,mf1y,mf1z,msize,'mutual friction')
%return
mesh_spectrum(mf1x,mf1y,mf1z,msize,'mutual friction',0)