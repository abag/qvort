%function andrew2pov(filenumber)
%Converts the vortex data used by Andrew Baggaley to a PovRay
%format into the object "Vortices" in a file "mesh.pov" that can be 
%loaded to make 3D figures with including the following lines:
%
%Also plots the configuration with simple plot3 command in matlab.  
%Run this with a filenumber as an input, e.g. andrew2pov(100)
%NOTE: Currently scales the vortex coordinates by 100!
function andrew2pov(filenumber)
%check filenumber has been set
if exist('filenumber')==0
  disp('you have not set filenumber')
  disp('aborting code and type "help vortex_plot" for more options')
  return
end
filename=sprintf('data/var%04d.log',filenumber);
%we set the dimensions of the box here
%this is overridden if we have periodic B.C.'s
box_size=.005 ;

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
  u2=fread(fid,number_of_particles,'float64');
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
    u2(j)=dummy_vect(6);
  end
  f=uint16(f);
end

%now create vectors to plot%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clf
for j=1:number_of_particles
  if round(f(j))==0
  else
    dummy_x(1,1)=x(j);
    dummy_x(2,1)=x(round(f(j)));
    dummy_x(1,2)=y(j);
    dummy_x(2,2)=y(round(f(j)));
    dummy_x(1,3)=z(j);
    dummy_x(2,3)=z(round(f(j)));
    dist=sqrt((dummy_x(1,1)-dummy_x(2,1))^2+(dummy_x(1,2)-dummy_x(2,2))^2+(dummy_x(1,3)-dummy_x(2,3))^2);
    if (dist<0.5*dims(2))
      plot3(dummy_x(1:2,1),dummy_x(1:2,2),dummy_x(1:2,3),'b-')
      hold on
    end %if (dist<0.5*dims(2))
  end %if round(f(j))==0
end

%The following contains the 

%now create the povray script %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%wmax=max(W); %maximum value for the smoothed vorticity
tones = 1000; %number of different color tones
ColMap=colormap(jet(tones+1));
%tuberad=0.02;
%tuberad=0.02/10;
tuberad=0.000025;
scale = 1000.00; %Added in order to get povray working correctly!

tuberad = scale*tuberad;
fid = fopen('mesh.pov','w');
fprintf(fid,'#declare Vortices = \n');
fprintf(fid,'union {\n');
for j=1:number_of_particles
  if round(f(j))==0
  else
    dummy_x(1,1)=x(j);
    dummy_x(2,1)=x(round(f(j)));
    dummy_x(1,2)=y(j);
    dummy_x(2,2)=y(round(f(j)));
    dummy_x(1,3)=z(j);
    dummy_x(2,3)=z(round(f(j)));
    dist=sqrt((dummy_x(1,1)-dummy_x(2,1))^2+(dummy_x(1,2)-dummy_x(2,2))^2+(dummy_x(1,3)-dummy_x(2,3))^2);
    if (dist<0.5*dims(2))
      %wval=W(j); %use when knowing the smoothed circulation
      %if (onlycoherent == 0) | (wval>treshold*wrms) 
        fprintf(fid,'sphere_sweep {\n');
        fprintf(fid,'linear_spline\n');
        N = 2; %linear interpolation
        fprintf(fid,'%5i\n',N);
        dummy_x = scale*dummy_x;
        fprintf(fid,'<%10.4f, %10.4f, %10.4f>, %10.4f\n',[dummy_x(1,1) dummy_x(1,2) dummy_x(1,3) tuberad]);  
        fprintf(fid,'<%10.4f, %10.4f, %10.4f>, %10.4f\n',[dummy_x(2,1) dummy_x(2,2) dummy_x(2,3) tuberad]);
        %fprintf(fid,'tolerance 0.1\n');
        %colorInd=round(tones*wval/wmax)+1;    
        %rgbvals=ColMap(colorInd,:);
        fprintf(fid,'pigment { color red 0 green 1 blue 0 }\n');
        %fprintf(fid,'pigment { color red %6.4f green %6.4f blue %6.4f }\n',rgbvals);
        fprintf(fid,'finish { ambient 0.2 diffuse 0.99 phong 1 }\n');
        fprintf(fid,'}\n');
      %end %if
      
    end %if (dist<0.5*dims(2))
  end %if round(f(j))==0
end
fprintf(fid,'}\n');%end of union




%Next make the colorbar (cylinder of height 1 and radius 0.03)
%using the colormap defined above.
fprintf(fid,'#declare ColorBar = \n');
fprintf(fid,'union {\n');
for i = 1:size(ColMap,1)
  rgbvals=ColMap(i,:);
  fprintf(fid,'cylinder {\n');
  z1=(i-1)/(tones+1)-0.5;
  z2=i/(tones+1)-0.5;
  rr1 = [0 0 z1];
  rr2 = [0 0 z2];
  rt  = 0.03;
  fprintf(fid,'<%10.4f, %10.4f, %10.4f>, <%10.4f, %10.4f, %10.4f>, %10.4f\n',[rr1 rr2 rt]);  
%  fprintf(fid,'sphere_sweep {\n');
%  fprintf(fid,'linear_spline\n');
%  N = 2;
%  fprintf(fid,'%5i\n',N);  
%  fprintf(fid,'<%10.4f, %10.4f, %10.4f>, %10.4f\n',[rr1 rt]);  
%  fprintf(fid,'<%10.4f, %10.4f, %10.4f>, %10.4f\n',[rr2 rt]);
  fprintf(fid,'pigment { color red %6.4f green %6.4f blue %6.4f }\n',rgbvals);
  fprintf(fid,'finish { ambient 0.2 diffuse 0.99 phong 1 }\n');
  fprintf(fid,'}\n');
end
fprintf(fid,'}\n');%end of union
  


