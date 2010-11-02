%not ready yet
function corrdim=corr_dim(filenumber,option)
if nargin==1     
  option='plot';
end
filename=sprintf('data/var%04d.log',filenumber);
%set options based on varargin
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%get the dimensions information from dims.log
dims=load('./data/dims.log');
delta=dims(1);  %set the resolution
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
%now we need to loop over particles and calculate correlation dimension
dist(1:number_of_particles,1:number_of_particles)=100.;
for i=1:number_of_particles
  if f(i)==0  ; continue ; end
  for j=1:number_of_particles
    if f(j)==0 ; continue ; end
    if i==j ; continue ; end
    next=i ;
    same_loop=0;
    for k=1:number_of_particles
      next=f(next);  
      if next==j ; same_loop=1 ; end
      if next==i ; break ; end
    end
    if  same_loop==0 ; continue ; end
    dist(i,j)=sqrt((x(i)-x(j))^2+(y(i)-y(j))^2+(z(i)-z(j))^2);
  end
end
switch option
  case 'plot'
    figure('Name','distance array')
    pcolor(dist/delta) ; shading interp ; colorbar
    set(gca,'FontSize',14);
    xlabel('pcount index','FontSize',12);
    ylabel('pcount index','FontSize',12);
    pause
end
for i=1:15
    rad(i)=delta*i;
    corr(i)=sum(sum(dist<rad(i)));
end
switch option
  case 'plot'
    plot(rad,corr,'LineWidth',2)
    set(gca,'FontSize',14);
    xlabel('radius','FontSize',12);
    ylabel('particle count','FontSize',12);
end
p=polyfit(log(rad),log(corr),1);
corrdim=p(1);
end
