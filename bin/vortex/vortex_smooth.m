function slope=vortex_smooth(filenumber,skip,option)
if nargin==2     
  option='plot';
end
filename=sprintf('data/var%04d.log',filenumber);
%set options based on varargin
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
%now order these particles, start with i=1
newx(1,1)=x(200);
newx(2,1)=y(200);
newx(3,1)=z(200);
next=f(200);
counter=2;
for i=2:number_of_particles
    newx(1,counter)=x(next);   
    newx(2,counter)=y(next);    
    newx(3,counter)=z(next);
    next=f(next);
    if (next==200)
        break
    end
    counter=counter+1;
end
s=size(newx);
counter=1;
for i=1:skip:s(2)
    newx2(1,counter)=newx(1,i);
    newx2(2,counter)=newx(2,i);
    newx2(3,counter)=newx(3,i);
    counter=counter+1;
end
%plot3(newx2(1,:),newx2(2,:),newx2(3,:),'b-',newx(1,:),newx(2,:),newx(3,:),'ro');
%pause
%now we want to use cubic spline
s2=size(newx2);
F=spline((1:s2(2)),newx2);
step=s2(2)/s(2);
step=0.05;
t=[1:step:s2(2)];
Ft=ppval(F,t);
%now find distances
switch option
  case 'plot'
     plot3(Ft(1,:),Ft(2,:),Ft(3,:),'b-','LineWidth',1.5)
     hold on
     plot3(newx(1,:),newx(2,:),newx(3,:),'ro');
     hold off
     set(gca,'FontSize',14)
     xlabel('x','FontSize',14)
     ylabel('y','FontSize',14)
     zlabel('z','FontSize',14)
     pause
end
dist(1:(s(2)-10))=10.;
s2=size(Ft);
len(1:(s(2)-10))=0;
for i=1:(s(2)-10)
    for j=1:s2(2)
        dum_len=sqrt((newx(1,i)-newx(1,i+1))^2+(newx(2,i)-newx(2,i+1))^2+(newx(3,i)-newx(3,i+1))^2);
        dum_dist=sqrt((newx(1,i)-Ft(1,j))^2+(newx(2,i)-Ft(2,j))^2+(newx(3,i)-Ft(3,j))^2);
        if dum_dist<dist(i)
          dist(i)=dum_dist;
          if i==1
            len(i)=dum_len;
          else
            len(i)=len(i-1)+dum_len;
          end
        end
    end
end

%plot(len,dist)
%pause
p=abs(fft(dist))/((s(2)-10)/2) ;
p=p(1:floor((s(2)-10)/2)).^2;
k=(1:floor((s(2)-10)/2))/max(len);
k2=1.:0.5:162;
p2=interp1(k,p,k2);
dmin=floor((s(2)-10)/3);
dmax=floor((s(2)-10)/2);
dump=p(dmin:dmax);
dumk=k(dmin:dmax);
s=size(k);
slope=polyfit(log(k(60:s(2))),log(p(60:s(2))),1);
dummp=log(k)*slope(1)+slope(2);
switch option
  case 'plot'
    plot(log(k),log(p),'k','LineWidth',2) ; hold on
    plot(log(k),dummp,'--r','LineWidth',1.5); hold off
     set(gca,'FontSize',14)
     xlabel('log k','FontSize',14)
     ylabel('log A(k)','FontSize',14)
end
end
