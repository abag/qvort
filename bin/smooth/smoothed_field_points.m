function mesh_plot(filenumber,filenumber2,firecolor)
if nargin<3
  firecolor=0;
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
fclose(fid)
zslice(1:msize,1:msize)=wz(msize/2,:,:);
imagesc(interp(x,2),interp(x,2),interp2(zslice,2));
if firecolor==1
  colormap(fireprint)
end
colorbar
hold on
%now load in vortex points
filename=sprintf('data/var%04d.log',filenumber2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%get the dimensions information from dims.log
dims=load('./data/dims.log');
if dims(4)==1
  fid=fopen(filename);
  if fid<0
      disp('var file does not exist, exiting script')
      return
  end
  time=fread(fid,1,'float64');
  number_of_particles=fread(fid,1,'int');
  x=fread(fid,number_of_particles,'float64');
  y=fread(fid,number_of_particles,'float64');
  z=fread(fid,number_of_particles,'float64');
  f=fread(fid,number_of_particles,'int');
  u=fread(fid,number_of_particles,'float64');
else
  fid=fopen(filename);
  if fid<0
      disp('var file does not exist, exiting script')
      return
  end
  %read the time
  tline=fgetl(fid);
  dummy=textscan(tline, '%f');
  time=dummy{:};
  %how many particles
  tline=fgetl(fid);
  dummy=textscan(tline, '%d');
  number_of_particles=dummy{:};
  %get the particles%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  for j=1:number_of_particles
    tline=fgetl(fid);
    dummy=textscan(tline, '%f');
    dummy_vect=dummy{:};
    x(j)=dummy_vect(1);
    y(j)=dummy_vect(2);
    z(j)=dummy_vect(3);
    f(j)=dummy_vect(4);
    u(j)=dummy_vect(5);
  end
  f=uint16(f);
end
%now create vectors to plot%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j=1:number_of_particles
  if round(f(j))==0
  else
    if z(j)*z(round(f(j)))<0
      testd=abs(z(j)-z(round(f(j))));
      if testd<dims(2)/10
        if z(j)>z(round(f(j)))
          plot(x(j),y(j),'k*')
        else
          plot(x(j),y(j),'ko')
        end 
        hold on
      end
    end
  end
end
xlabel('x','FontSize',16)
ylabel('y','FontSize',16)
set(gca,'FontSize',16)
