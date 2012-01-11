function velocity_snapshot(filenumber)
  filename=sprintf('data/uu%04d.dat',filenumber);
  fid=fopen(filename);
  if fid<0
    disp('var file does not exist, exiting script')
    return
  end
  %disp(sprintf('reading var %04d',i))
  time=fread(fid,1,'float64');
  number_of_particles=fread(fid,1,'int');
  ux=fread(fid,number_of_particles,'float64');
  uy=fread(fid,number_of_particles,'float64');
  uz=fread(fid,number_of_particles,'float64');
  subplot(1,3,1)
    ksdensity(ux)
  subplot(1,3,2)
    ksdensity(uy) 
  subplot(1,3,3)
    ksdensity(uz)
end
