function vortex_point_slice(filenumber)
%load in vortex points
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
%now create vectors to plot%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
counterdown=1;
counterup=1;
for i=1:500
plane=-dims(2)/3. + (2*dims(2)/3)*rand;
for j=1:number_of_particles
  if round(f(j))==0
  else
    testd=sqrt((x(j)-x(round(f(j))))^2+(y(j)-y(round(f(j))))^2+(z(j)-z(round(f(j))))^2);
    if (testd<10.*dims(1))
          if (x(j)<plane) && (x(round(f(j)))>plane)
              pointup(counterup,1)=y(j);
              pointup(counterup,2)=z(j);
              counterup=counterup+1;
          elseif (x(j)>plane) && (x(round(f(j)))<plane)
              pointdown(counterdown,1)=y(j);
              pointdown(counterdown,2)=z(j);
              counterdown=counterdown+1;
          end
    end
  end
end
[theta1 r1]=cart2pol(pointup(:,1),pointup(:,2));
[theta2 r2]=cart2pol(pointdown(:,1),pointdown(:,2));
ind1=find(r1<dims(2)/2);
ind2=find(r2<dims(2)/2);
Gamma_avg(i)=(size(r1(ind1),1)-size(r2(ind2),1))^2;
%if Gamma_avg(i)>36
%polar(theta1(ind1),r1(ind1),'o')
%hold on
%polar(theta2(ind2),r2(ind2),'+')
%hold off
%pause
%end
clear r1 r2 theta1 theta2 pointup pointdown ind1 ind2
end
mean(Gamma_avg)
max(Gamma_avg)
ts=load('./data/ts.log');
ell=1./sqrt(ts(1,6)/dims(2)^3);
log(mean(Gamma_avg))/log(dims(2)/ell)
