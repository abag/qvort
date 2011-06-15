function [f xi]=vortex_dist_hist(filenumber,option)
if nargin==1
  option='plot';
end
filename=sprintf('data/var%04d.log',filenumber);
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
fclose(fid);
number_of_particles=uint16(number_of_particles);
sep(number_of_particles*number_of_particles/2)=0.;
%now create vectors to plot%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
counter=1;
for j=1:number_of_particles
  if round(f(j))==0
  else
    for i=j:number_of_particles
      if round(f(i))==0
      else
        sep(counter)=sqrt((x(j)-x(i))^2+(y(j)-y(i))^2+(z(j)-z(i))^2);
        counter=counter+1;
      end
    end
  end
end
disp('finished computing distances')
[f xi]=ksdensity(sep,'support',positive);
switch option
  case 'plot'
   figure('Name','Histogram of particle separation')
     hist(sep);
     xlabel('separation','FontSize',16);
     ylabel('N','FontSize',16);
     set(gca,'FontSize',16)
   figure('Name','PDF particle separation') 
     plot(xi,f,'-k','LineWidth',2)
     xlabel('separation','FontSize',16)
     ylabel('PDF','FontSize',16)
     set(gca,'FontSize',16)
end
clear sep