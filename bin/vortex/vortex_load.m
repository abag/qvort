function vortex_load(filenumber)
global dims box_size
global x y z
global f u u_mf v_curv v_stretch f_mf
global number_of_particles
global u_mf_x u_mf_y u_mf_z
global ux uy uz
global v_f_mf_x v_f_mf_y v_f_mf_z
%check filenumber has been set
if exist('filenumber')==0
  disp('you have not set filnumber')
  disp('aborting code and type "help vortex_plot" for more options')
  return
end
filename=sprintf('data/var%04d.log',filenumber);
%get the dimensions information from dims.log
dims=load('./data/dims.log');
box_size=dims(2);
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
ux=fread(fid,number_of_particles,'float64');
uy=fread(fid,number_of_particles,'float64');
uz=fread(fid,number_of_particles,'float64');
u_mf_x=fread(fid,number_of_particles,'float64');
u_mf_y=fread(fid,number_of_particles,'float64');
u_mf_z=fread(fid,number_of_particles,'float64');
v_curv=fread(fid,number_of_particles,'float64');
v_stretch=fread(fid,number_of_particles,'float64');
v_f_mf_x=fread(fid,number_of_particles,'float64');
v_f_mf_y=fread(fid,number_of_particles,'float64');
v_f_mf_z=fread(fid,number_of_particles,'float64');
fclose(fid);
%now create arrays for u and u_mf
u=sqrt(ux.^2+uy.^2+uz.^2);
u_mf=sqrt(u_mf_x.^2+u_mf_y.^2+u_mf_z.^2);
f_mf=sqrt(v_f_mf_x.^2+v_f_mf_y.^2+v_f_mf_z.^2);



