function [t ld_inside]=plot_vortex_within_smooth(filenumber,filenumber2)
filename=sprintf('data/smoothed_field%03d.dat',filenumber);
load ./data/sm_dims.log;
msize=sm_dims(1);
fid=fopen(filename);
if fid<0
  disp('file does not exist, exiting script')
  return
end
t=fread(fid,1,'float64');
xmesh=fread(fid,msize,'float64');
wx=fread(fid,msize^3,'float64');
wy=fread(fid,msize^3,'float64');
wz=fread(fid,msize^3,'float64');
rms_w=sqrt(mean(wx.^2+wy.^2+wz.^2))
wx=reshape(wx,msize,msize,msize);
wy=reshape(wy,msize,msize,msize);
wz=reshape(wz,msize,msize,msize);
mod_w=sqrt(wx.^2+wy.^2+wz.^2);
%mod_w=permute(mod_w,[2 1 3]);
fclose(fid);
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
[sw3,sw2,sw1]=ind2sub(size(mod_w),find(mod_w>2.2*rms_w));
total_volume=length(sw1)*((dims(2)/msize)^3)
total_length=0. ;
for j=1:number_of_particles
  if round(f(j))==0
  else
    dummy_x(1,1)=x(j);
    dummy_x(2,1)=x(round(f(j)));
    dummy_x(1,2)=y(j);
    dummy_x(2,2)=y(round(f(j)));
    dummy_x(1,3)=z(j);
    dummy_x(2,3)=z(round(f(j)));
    can_plot=0;
    for i=1:length(sw1)
      if sw1(i)==msize || sw2(i)==msize || sw3(i)==msize
        continue
      end
      if sw1(i)==msize-1 || sw2(i)==msize-1 || sw3(i)==msize-1
        continue
      end
      if sw1(i)==1 || sw2(i)==1 || sw3(i)==1
        continue
      end
      if dummy_x(1,1)>xmesh(sw1(i)-1) && dummy_x(1,1)<xmesh(sw1(i)+2) &&        dummy_x(1,2)>xmesh(sw2(i)-1) && dummy_x(1,2)<xmesh(sw2(i)+2) && dummy_x(1,      3)>xmesh(sw3(i)-1) && dummy_x(1,3)<xmesh(sw3(i)+2)
        can_plot=1;
        break
      end
    end
    if can_plot==1
      dist=sqrt((dummy_x(1,1)-dummy_x(2,1))^2+(dummy_x(1,2)-dummy_x(2,2))^2+(dummy_x(1,3)-dummy_x(2,3))^2);
      if dist<dims(2)/4
        total_length=total_length+dist;
      end
    end 
  end
end
ld_inside=total_length/total_volume;
