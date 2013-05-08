function vortex_load(filenumber)
global dims
global x y z
global f
global number_of_particles
%check filenumber has been set
if exist('filenumber')==0
  disp('you have not set filnumber')
  disp('aborting code and type "help vortex_plot" for more options')
  return
end
filename=sprintf('data/var_interp%04d.log',filenumber);
%get the dimensions information from dims.log
dims=load('./data/dims.log');
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
