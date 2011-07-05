function run_comp_all(filenumbers)
for i=filenumbers
filename1=sprintf('./run1/data/var%04d.log',i);
dims=load('./run1/data/dims.log');
if dims(4)==1
  fid=fopen(filename1);
  if fid<0
      disp('var file does not exist, exiting script')
      return
  end
  time=fread(fid,1,'float64');
  number_of_particles=fread(fid,1,'int');
  x1=fread(fid,number_of_particles,'float64');
  y1=fread(fid,number_of_particles,'float64');
  z1=fread(fid,number_of_particles,'float64');
  f1=fread(fid,number_of_particles,'int');
  u1=fread(fid,number_of_particles,'float64');
  u21=fread(fid,number_of_particles,'float64');
else 
  fid=fopen(filename1);
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
    x1(j)=dummy_vect(1);
    y1(j)=dummy_vect(2);
    z1(j)=dummy_vect(3);
    f1(j)=dummy_vect(4);
    u1(j)=dummy_vect(5);
    u21(j)=dummy_vect(6);
  end
  f1=uint16(f1);
end
filename2=sprintf('./run2/data/var%04d.log',i);
dims=load('./run2/data/dims.log');
if dims(4)==1
  fid=fopen(filename2);
  if fid<0
      disp('var file does not exist, exiting script')
      return
  end
  true_time(i)=fread(fid,1,'float64');
  number_of_particles=fread(fid,1,'int');
  x2=fread(fid,number_of_particles,'float64');
  y2=fread(fid,number_of_particles,'float64');
  z2=fread(fid,number_of_particles,'float64');
  f2=fread(fid,number_of_particles,'int');
  u2=fread(fid,number_of_particles,'float64');
  u22=fread(fid,number_of_particles,'float64');
else 
  fid=fopen(filename2);
  if fid<0
      disp('var file does not exist, exiting script')
      return
  end
  %read the time
  tline=fgetl(fid);
  dummy=textscan(tline, '%f');
  true_time(i)=dummy{:};
  %how many particles
  tline=fgetl(fid);
  dummy=textscan(tline, '%d');  
  number_of_particles=dummy{:};
  %get the particles%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  for j=1:number_of_particles
    tline=fgetl(fid);
    dummy=textscan(tline, '%f');
    dummy_vect=dummy{:};
    x2(j)=dummy_vect(1);
    y2(j)=dummy_vect(2);
    z2(j)=dummy_vect(3);
    f2(j)=dummy_vect(4);
    u2(j)=dummy_vect(5);
    u22(j)=dummy_vect(6);
  end
  f2=uint16(f2);
end
filename3=sprintf('./run3/data/var%04d.log',i);
dims=load('./run3/data/dims.log');
if dims(4)==1
  fid=fopen(filename3);
  if fid<0
      disp('var file does not exist, exiting script')
      return
  end
  time=fread(fid,1,'float64');
  number_of_particles=fread(fid,1,'int');
  x3=fread(fid,number_of_particles,'float64');
  y3=fread(fid,number_of_particles,'float64');
  z3=fread(fid,number_of_particles,'float64');
  f3=fread(fid,number_of_particles,'int');
  u3=fread(fid,number_of_particles,'float64');
  u32=fread(fid,number_of_particles,'float64');
else 
  fid=fopen(filename3);
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
    x3(j)=dummy_vect(1);
    y3(j)=dummy_vect(2);
    z3(j)=dummy_vect(3);
    f3(j)=dummy_vect(4);
    u3(j)=dummy_vect(5);
    u32(j)=dummy_vect(6);
  end
  f3=uint16(f3);
end
diff(i)=mean(sqrt((x1-x2).^2+(y1-y2).^2+(z1-z2).^2));
diff2(i)=mean(sqrt((x1-x3).^2+(y1-y3).^2+(z1-z3).^2));
end
plot(true_time,diff/dims(1),true_time,diff2/dims(1))