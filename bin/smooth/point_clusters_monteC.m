function point_clusters_monteC(filenumber,iteration)
%check nargin
if nargin<2
  iteration=1;
end
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
counter1(1:iteration)=1; %set counter to 1
counter2(1:iteration)=1; %set counter to 1
main_counter=1
locs=linspace(dims(1),dims(2)/3,30);
K_mean(1:2*iteration,1:30)=0.;
Kmax=0. ; 
for i=1:iteration
  %pick a random plane in xy
  rand_plane=rand*dims(2)-dims(2)/2.;
  for j=1:number_of_particles
    if round(f(j))==0
      %skip empty particles
    else
      if z(j)<rand_plane && z(round(f(j)))>rand_plane
        field(2*i-1,counter1(i),1)=x(j);
        field(2*i-1,counter1(i),2)=y(j);
        counter1(i)=counter1(i)+1;
      end
      if z(j)>rand_plane && z(round(f(j)))<rand_plane
        field(2*i,counter2(i),1)=x(j);
        field(2*i,counter2(i),2)=y(j);
        counter2(i)=counter2(i)+1;
      end
    end
  end
%   subplot(2,1,1)
%   plot(field(i,:,1),field(i,:,2),'o')
%   subplot(2,1,2)
%   plot(field(i+1,:,1),field(i+1,:,2),'o')
%   pause
  K = ripleyK(squeeze(field(i,:,1:2)),locs,[-dims(2)/2. dims(2)/2. -dims(2)/2. dims(2)/2. -dims(2)/2. dims(2)/2.]);
  if max(sqrt(K/pi)-locs')<0.02
    plot(locs,sqrt(K/pi)-locs') ; hold on
    K_mean(main_counter,:)=sqrt(K/pi)-locs';
    main_counter=main_counter+1;
  end
  K = ripleyK(squeeze(field(i+1,:,1:2)),locs,[-dims(2)/2. dims(2)/2. -dims(2)/2. dims(2)/2. -dims(2)/2. dims(2)/2.]);
  if max(sqrt(K/pi)-locs')<0.02
    plot(locs,sqrt(K/pi)-locs') ; hold on  
    K_mean(main_counter,:)=sqrt(K/pi)-locs';
    main_counter=main_counter+1;
  end
  if sum(K)>Kmax
    iKmax=i;
    Kmax=sum(K);
  end
end
figure
errorbar(locs,mean(K_mean),std(K_mean))
%%%%%%%%%%%%%%%%%%%%%%%%%%RANDOM DATA%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear K K_mean
for i=1:200
  dumx(:,1)=rand(300,1)*dims(2)-dims(2)/2.;
  dumx(:,2)=rand(300,1)*dims(2)-dims(2)/2.;
  K = ripleyK(dumx,locs,[-dims(2)/2. dims(2)/2. -dims(2)/2. dims(2)/2. -dims(2)/2. dims(2)/2.]);
  K_mean(i+1,:)=sqrt(K/pi)-locs';  
end
hold on
errorbar(locs,mean(K_mean),2*std(K_mean))
%figure
%  subplot(2,1,1)
%  plot(field(iKmax,:,1),field(iKmax,:,2),'o')
%  subplot(2,1,2)
%  plot(field(iKmax+1,:,1),field(iKmax+1,:,2),'o')
