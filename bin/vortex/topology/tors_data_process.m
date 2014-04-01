function tors_data_process(start,finish)
global tors
if nargin<2
    finish=start;
end
if nargin<1
    start=1;
    finish=start;
end
for i=start:finish
  filename=sprintf('data/torsion_data%04d.dat',i);
  fid=fopen(filename);
  if fid<0
    disp('var file does not exist, exiting script')
    return
  end
  %disp(sprintf('reading var %04d',i))
  time=fread(fid,1,'float64');
  number_of_particles=fread(fid,1,'int');
  if i==start
    tors=fread(fid,number_of_particles,'float64');
  else
    l=length(tors);
    tors(l+1:l+number_of_particles)=fread(fid,number_of_particles,'float64');
  end
  fclose(fid);
end
