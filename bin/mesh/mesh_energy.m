function mesh_energy(start,finish,skip)
if nargin<1
  disp('I at least need finish filenumbers')
  return
elseif nargin<2
  disp('Assuming number given is final, start set to 1')
  final=start;
  start=1;
  skip=1.
elseif nargin<3
  disp('skip set to 1')
  skip=1;
end
load data/dims.log;
bsize=dims(2); %box size
msize=dims(3); %mesh size
if (msize==0) 
  disp('mesh size is zero exiting script')
  return
end
for i=start:skip:finish
  i
  filename=sprintf('data/mesh%03d.dat',i);
  fid=fopen(filename);
  if fid<0
    disp('mesh file does not exist, exiting script')
    return
  end
  t(i)=fread(fid,1,'float64');
  x=fread(fid,msize,'float64');
  unormx=fread(fid,msize^3,'float64');
  unormy=fread(fid,msize^3,'float64');
  unormz=fread(fid,msize^3,'float64');
  ux=fread(fid,msize^3,'float64');
  uy=fread(fid,msize^3,'float64');
  uz=fread(fid,msize^3,'float64');
  %%%%%%%%%%%%%CHECKING FOR ABNORMAL VELOCITY%%%%%%%%%%%%%%%%
  uu=sqrt(ux.^2+uy.^2+uz.^2);
  vcoff=50. ;
  index = find(uu > vcoff);
  ux(index)=0. ; uy(index)=0. ;uz(index)=0. ;
  clear index
  u2(i)=sum(sum(sum(ux.^2+uy.^2+uz.^2)))/(2.*(msize/bsize)^3);
  fclose(fid);
end
disp('saving t and u2 to mesh_energy.mat')
save mesh_energy.mat t u2
disp('use mesh_energy_fit to plot')
