function vortex_LxLyLz
%load in time series data to get number of files
A=load('./data/ts.log');
%get the dimensions information from dims.log
dims=load('./data/dims.log');
%now load in vortex points
for i=1:size(A,1)
filename=sprintf('./data/var%04d.log',i);
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
  fclose(fid);
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
  fclose(fid);
end
total_Lx(i)=0. ;
total_Ly(i)=0. ;
total_Lz(i)=0. ;
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
    dist=sqrt((dummy_x(1,1)-dummy_x(2,1))^2+(dummy_x(1,2)-dummy_x(2,2))^2+(dummy_x(1,3)-dummy_x(2,3))^2);
    if (dist<dims(2)/4.)
      total_Lx(i)=total_Lx(i)+abs(dummy_x(1,1)-dummy_x(2,1));
      total_Ly(i)=total_Ly(i)+abs(dummy_x(1,2)-dummy_x(2,2));
      total_Lz(i)=total_Lz(i)+abs(dummy_x(1,3)-dummy_x(2,3));
    end 
  end
end
end
subplot(3,1,1)
  plot(A(:,2),total_Lx,'k','LineWidth',2)
  xlabel('t','FontSize',16)
  ylabel('L_x','FontSize',16)
  set(gca,'FontSize',16)
subplot(3,1,2)
  plot(A(:,2),total_Ly,'r','LineWidth',2)
  xlabel('t','FontSize',16)
  ylabel('L_y','FontSize',16)
  set(gca,'FontSize',16)
subplot(3,1,3)
  plot(A(:,2),total_Lz,'b','LineWidth',2)
  xlabel('t','FontSize',16)
  ylabel('L_z','FontSize',16)
  set(gca,'FontSize',16)
  
 t_L=A(:,2);
save ./LxLyLz.mat t_L total_Lx total_Ly total_Lz
