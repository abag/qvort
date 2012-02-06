function vortex_point_slice(filenumber,plane)
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
counterx=1;
countery=1;
counterz=1; 
for j=1:number_of_particles
  if round(f(j))==0
  else
    testd=sqrt((x(j)-x(round(f(j))))^2+(y(j)-y(round(f(j))))^2+(z(j)-z(round(f(j))))^2);
    if (testd<8.*dims(1))
      switch plane
        case 'x'
          if x(j)*x(round(f(j)))<0
              pointx(counterx,1)=y(j);
              pointx(counterx,2)=z(j);
              counterx=counterx+1;
          end
        case 'y'
          if y(j)*y(round(f(j)))<0
              pointy(countery,1)=z(j);
              pointy(countery,2)=x(j);
              countery=countery+1;
          end
        case 'z'
          if z(j)*z(round(f(j)))<0
              pointz(counterz,1)=x(j);
              pointz(counterz,2)=y(j);
              counterz=counterz+1;
          end
        case 'yz'
          if y(j)*y(round(f(j)))<0
              pointy(countery,1)=z(j);
              pointy(countery,2)=x(j);
              countery=countery+1;
          end
          if z(j)*z(round(f(j)))<0
              pointz(counterz,1)=x(j);
              pointz(counterz,2)=y(j);
              counterz=counterz+1;
          end
      end
    end
  end
end
switch plane
  case 'x'
    plot(pointx(:,1),pointx(:,2),'ko','MarkerFaceColor','k')
    axis([-dims(2)/2 dims(2)/2 -dims(2)/2 dims(2)/2])
    set(gca,'FontSize',16)
    xlabel('y','FontSize',14)
    ylabel('z','FontSize',14)
  case 'y'
    plot(pointy(:,1),pointy(:,2),'k.')
    axis([-dims(2)/2 dims(2)/2 -dims(2)/2 dims(2)/2])
    set(gca,'FontSize',16)
    xlabel('z','FontSize',14)
    ylabel('x','FontSize',14)
  case 'z'
    plot(pointz(:,1),pointz(:,2),'k.')
    axis([-dims(2)/2 dims(2)/2 -dims(2)/2 dims(2)/2])
    set(gca,'FontSize',16)
    xlabel('x','FontSize',14)
    ylabel('y','FontSize',14)
  case 'yz'
	subplot(3,4,[1 5 9])
    	plot(pointy(:,1),pointy(:,2),'k.')
        axis([-.2 .2 -dims(2)/2 dims(2)/2],'equal')
        set(gca,'FontSize',16)
    	xlabel('z','FontSize',16)
    	ylabel('x','FontSize',16)
	subplot(3,4,[6 7 8])
    	plot(pointz(:,1),pointz(:,2),'k.')
        axis([-dims(2)/2 dims(2)/2 -.2 .2],'equal')
        set(gca,'FontSize',16) 
    	xlabel('x','FontSize',16)
    	ylabel('y','FontSize',16)
end
s1='t=';
s2=num2str(time);
str=strcat(s1,s2);
if (dims(2)>0.)
  text(-.6*dims(2),.6*dims(2),.7*dims(2),str,'FontSize',16)
else
  text(-0.06,0.06,0.07,str,'FontSize',16)
end
