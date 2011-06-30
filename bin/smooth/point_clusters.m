function point_clusters(filenumber)
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
counterx1=1; counterx2=1;
countery1=1; countery2=1;
counterz1=1; counterz2=1;
for j=1:number_of_particles
  if round(f(j))==0
  else
    testd=sqrt((x(j)-x(round(f(j))))^2+(y(j)-y(round(f(j))))^2+(z(j)-z(round(f(j))))^2);
    if (testd<8.*dims(1))
      if x(j)*x(round(f(j)))<0
        if x(j)>x(round(f(j)))
          pointupx(counterx1,1)=z(j);
          pointupx(counterx1,2)=y(j);
          counterx1=counterx1+1;
        else
          pointdownx(counterx2,1)=z(j);
          pointdownx(counterx2,2)=y(j);
          counterx2=counterx2+1;
        end 
      end
      if y(j)*y(round(f(j)))<0
        if y(j)>y(round(f(j)))
          pointupy(countery1,1)=x(j);
          pointupy(countery1,2)=z(j);
          countery1=countery1+1;
        else
          pointdowny(countery2,1)=x(j);
          pointdowny(countery2,2)=z(j);
          countery2=countery2+1;
        end 
      end
      if z(j)*z(round(f(j)))<0
        if z(j)>z(round(f(j)))
          pointupz(counterz1,1)=x(j);
          pointupz(counterz1,2)=y(j);
          counterz1=counterz1+1;
        else
          pointdownz(counterz2,1)=x(j);
          pointdownz(counterz2,2)=y(j);
          counterz2=counterz2+1;
        end 
      end
    end
  end
end
subplot(3,2,1)
  plot(pointupx(:,1),pointupx(:,2),'o')
subplot(3,2,2)
  plot(pointdownx(:,1),pointdownx(:,2),'o')
subplot(3,2,3)
  plot(pointupy(:,1),pointupy(:,2),'o')
subplot(3,2,4)
  plot(pointdowny(:,1),pointdowny(:,2),'o')
subplot(3,2,5)
  plot(pointupz(:,1),pointupz(:,2),'o')
subplot(3,2,6)
  plot(pointdownz(:,1),pointdownz(:,2),'o')
%%%%%%%%%%%%%%%%%%%%%%%%%%Now do ripley K%%%%%%%%%%%%%%%%%%%
locs=linspace(dims(1),dims(2)/3,30);
Kx1 = ripleyK(pointupx,locs,[-dims(2)/2. dims(2)/2. -dims(2)/2. dims(2)/2. -dims(2)/2. dims(2)/2.]);
Kx2 = ripleyK(pointdownx,locs,[-dims(2)/2. dims(2)/2. -dims(2)/2. dims(2)/2. -dims(2)/2. dims(2)/2.]);
Ky1 = ripleyK(pointupy,locs,[-dims(2)/2. dims(2)/2. -dims(2)/2. dims(2)/2. -dims(2)/2. dims(2)/2.]);
Ky2 = ripleyK(pointdowny,locs,[-dims(2)/2. dims(2)/2. -dims(2)/2. dims(2)/2. -dims(2)/2. dims(2)/2.]);
Kz1 = ripleyK(pointupz,locs,[-dims(2)/2. dims(2)/2. -dims(2)/2. dims(2)/2. -dims(2)/2. dims(2)/2.]);
Kz2 = ripleyK(pointdownz,locs,[-dims(2)/2. dims(2)/2. -dims(2)/2. dims(2)/2. -dims(2)/2. dims(2)/2.]);
figure
plot(locs, sqrt(Kx1/pi)-locs',locs, sqrt(Kx2/pi)-locs',locs, sqrt(Ky1/pi)-locs',locs, sqrt(Ky2/pi)-locs',locs, sqrt(Kz1/pi)-locs',locs, sqrt(Kz2/pi)-locs')
X(:,1)=sqrt(Kx1/pi); X(:,2)=sqrt(Kx2/pi);
X(:,3)=sqrt(Ky1/pi); X(:,4)=sqrt(Ky2/pi);
X(:,5)=sqrt(Kz1/pi); X(:,6)=sqrt(Kz2/pi);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%SIMULATE RANDOM POINTS%%%%%%%%%%%%%%%%%%%%%
mulptiply=4;
pointupx(1:mulptiply*length(pointupx),:)=rand(mulptiply*length(pointupx),2)*dims(2)-dims(2)/2;
pointdownx(1:mulptiply*length(pointdownx),:)=rand(mulptiply*length(pointdownx),2)*dims(2)-dims(2)/2;
pointupy(1:mulptiply*length(pointupy),:)=rand(mulptiply*length(pointupy),2)*dims(2)-dims(2)/2;
pointdowny(1:mulptiply*length(pointdowny),:)=rand(mulptiply*length(pointdowny),2)*dims(2)-dims(2)/2;
pointupz(1:mulptiply*length(pointupz),:)=rand(mulptiply*length(pointupz),2)*dims(2)-dims(2)/2;
pointdownz(1:mulptiply*length(pointdownz),:)=rand(mulptiply*length(pointdownz),2)*dims(2)-dims(2)/2;
Kx1 = ripleyK(pointupx,locs,[-dims(2)/2. dims(2)/2. -dims(2)/2. dims(2)/2. -dims(2)/2. dims(2)/2.]);
Kx2 = ripleyK(pointdownx,locs,[-dims(2)/2. dims(2)/2. -dims(2)/2. dims(2)/2. -dims(2)/2. dims(2)/2.]);
Ky1 = ripleyK(pointupy,locs,[-dims(2)/2. dims(2)/2. -dims(2)/2. dims(2)/2. -dims(2)/2. dims(2)/2.]);
Ky2 = ripleyK(pointdowny,locs,[-dims(2)/2. dims(2)/2. -dims(2)/2. dims(2)/2. -dims(2)/2. dims(2)/2.]);
Kz1 = ripleyK(pointupz,locs,[-dims(2)/2. dims(2)/2. -dims(2)/2. dims(2)/2. -dims(2)/2. dims(2)/2.]);
Kz2 = ripleyK(pointdownz,locs,[-dims(2)/2. dims(2)/2. -dims(2)/2. dims(2)/2. -dims(2)/2. dims(2)/2.]);
Y(:,1)=sqrt(Kx1/pi); Y(:,2)=sqrt(Kx2/pi);
Y(:,3)=sqrt(Ky1/pi); Y(:,4)=sqrt(Ky2/pi);
Y(:,5)=sqrt(Kz1/pi); Y(:,6)=sqrt(Kz2/pi);
hold on
figure
errorbar(locs,locs-mean(Y'),std(Y'))
hold on
errorbar(locs,locs-mean(X'),std(X')/2)
%boundedline(locs,locs-mean(Y'),std(Y')','k-','alpha','transparency', 0.1)
%boundedline(locs,locs-mean(X'),std(X')/2,'k-','alpha')
