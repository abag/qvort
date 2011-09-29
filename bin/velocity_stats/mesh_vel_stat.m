function mesh_vel_stat(filenumbers,downsample,display)
if nargin<2
    downsample=0;
    display=0;
end
if nargin<3
    display=0;
else
    display=1;
    disp('producing iso surfaces before and after + paraview output')
end
load ./data/dims.log;
msize=dims(3);
if (msize==0) 
  disp('mesh size is zero exiting script')
  return
end
disp(sprintf('mesh size is: %04d',msize))
counter=0;
for i=filenumbers
  filename=sprintf('data/mesh%03d.dat',i);
  fid=fopen(filename);
  if fid<0
    disp('mesh file does not exist, exiting script')
    return
  end
  t=fread(fid,1,'float64');
  x=fread(fid,msize,'float64');
  dummy_unormx=fread(fid,msize^3,'float64');
  dummy_unormy=fread(fid,msize^3,'float64');
  dummy_unormz=fread(fid,msize^3,'float64');
  dummy_ux=fread(fid,msize^3,'float64');
  dummy_uy=fread(fid,msize^3,'float64');
  dummy_uz=fread(fid,msize^3,'float64');
  if downsample>0
      dummy_ux=reshape(dummy_ux,msize,msize,msize);
      dummy_uy=reshape(dummy_uy,msize,msize,msize);
      dummy_uz=reshape(dummy_uz,msize,msize,msize);
      if display==1
        mesh_iso(x,dummy_ux,dummy_uy,dummy_uz,msize,'super')
        writevtk(sqrt(dummy_ux.^2+dummy_uy.^2+dummy_uz.^2),'para_sup.vtk')
      end
      dummy_ux=smooth3(dummy_ux,'box',downsample);
      dummy_uy=smooth3(dummy_ux,'box',downsample);
      dummy_uz=smooth3(dummy_ux,'box',downsample);
      if display==1
        mesh_iso(x,dummy_ux,dummy_uy,dummy_uz,msize,'super_smooth')
        writevtk(sqrt(dummy_ux.^2+dummy_uy.^2+dummy_uz.^2),'para_sup_smooth.vtk')
      end
      dummy_ux=reshape(dummy_ux,msize^3,1);
      dummy_uy=reshape(dummy_uy,msize^3,1);
      dummy_uz=reshape(dummy_uz,msize^3,1);
  end
  ux(counter*msize^3+1:(counter+1)*msize^3)=dummy_ux;
  uy(counter*msize^3+1:(counter+1)*msize^3)=dummy_uy;
  uz(counter*msize^3+1:(counter+1)*msize^3)=dummy_uz;
  counter=counter+1;
end
save velocity.mat ux uy uz
