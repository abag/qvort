function KScompare(var1,var2)
filename1=sprintf('./run1/data/var%04d.log',var1);
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
filename2=sprintf('./run2/data/var%04d.log',var2);
dims=load('./run2/data/dims.log');
if dims(4)==1
  fid=fopen(filename2);
  if fid<0
      disp('var file does not exist, exiting script')
      return
  end
  time=fread(fid,1,'float64');
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
    x2(j)=dummy_vect(1);
    y2(j)=dummy_vect(2);
    z2(j)=dummy_vect(3);
    f2(j)=dummy_vect(4);
    u2(j)=dummy_vect(5);
    u22(j)=dummy_vect(6);
  end
  f2=uint16(f2);
end
diff=sqrt((x1-x2).^2+(y1-y2).^2+(z1-z2).^2)./sqrt(x1.^2+y1.^2+z1.^2);
dispme(1)=max(diff) ; dispme(2)=mean(diff) ; dispme(3)=min(diff) ;
dispme
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  filename1=sprintf('./run1/data/uu%04d.dat',var1);
  fid=fopen(filename1);
  if fid<0
    disp('var file does not exist, exiting script')
    return
  end
  time=fread(fid,1,'float64');
  number_of_particles=fread(fid,1,'int');
  ux1=fread(fid,number_of_particles,'float64');
  uy1=fread(fid,number_of_particles,'float64');
  uz1=fread(fid,number_of_particles,'float64');
  fclose(fid);
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  filename2=sprintf('./run2/data/uu%04d.dat',var2);
  fid=fopen(filename2);
  if fid<0
    disp('var file does not exist, exiting script')
    return
  end
  time=fread(fid,1,'float64');
  number_of_particles=fread(fid,1,'int');
  ux2=fread(fid,number_of_particles,'float64');
  uy2=fread(fid,number_of_particles,'float64');
  uz2=fread(fid,number_of_particles,'float64');
  fclose(fid);
  
diffu=sqrt((ux1-ux2).^2+(uy1-uy2).^2+(uz1-uz2).^2)./sqrt(ux1.^2+uy1.^2+uz1.^2);
dispme(1)=max(diffu) ; dispme(2)=mean(diffu) ; dispme(3)=min(diffu) ;
dispme
