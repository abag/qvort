function vel_stat_process(start,finish)
if nargin<2
    finish=start;
end
if nargin<1
    start=1;
    finish=start;
end
for i=start:finish
  filename=sprintf('data/uu%04d.dat',i);
  fid=fopen(filename);
  if fid<0
    disp('var file does not exist, exiting script')
    return
  end
  %disp(sprintf('reading var %04d',i))
  time=fread(fid,1,'float64');
  number_of_particles=fread(fid,1,'int');
  if i==start
    ux=fread(fid,number_of_particles,'float64');
    uy=fread(fid,number_of_particles,'float64');
    uz=fread(fid,number_of_particles,'float64');
  else
    l=length(ux);
    ux(l+1:l+number_of_particles)=fread(fid,number_of_particles,'float64');
    uy(l+1:l+number_of_particles)=fread(fid,number_of_particles,'float64');
    uz(l+1:l+number_of_particles)=fread(fid,number_of_particles,'float64');
  end
  fclose(fid);
end
markerx=1;
markery=1;
markerz=1;
u2=sqrt(ux.^2+uy.^2+uz.^2);
vcoff=5. ;
index = find(ux > vcoff);
ux(index) = [];
clear index;
index = find(uy > vcoff);
uy(index) = [];
clear index;
index = find(uz > vcoff);
uz(index) = [];
clear index ;
index = find(ux < -vcoff);
ux(index) = [];
clear index;
index = find(uy < -vcoff);
uy(index) = [];
clear index;
index = find(uz < -vcoff);
uz(index) = [];
clear index ;
index = find(u2 > vcoff);
u2(index) = [];
clear index;
ux=ux' ; uy=uy' ; uz=uz' ; u2=u2' ;
save velocity.mat ux uy uz u2
