function kdtree_povray_out(filenumber,plot,eval_rms)
global dims
global x y z
global f
global number_of_particles
global vort vort_rms
if nargin<2
  plot=1;
  eval_rms=2.;
elseif nargin<3
  eval_rms=2.;
end
%global test
vortex_load(filenumber);
counter=1;
for j=1:number_of_particles
  if f(j)==0 ; continue ; end
  dist=sqrt((x(j)-x(round(f(j))))^2+(y(j)-y(round(f(j))))^2+(z(j)-z(round(f(j))))^2);
  if (dist>0.5*dims(2)) ; continue ; end
  varray(counter,1)=x(j);
  varray(counter,2)=y(j);
  varray(counter,3)=z(j);
  varray(counter,4)=x(round(f(j)))-x(j);
  varray(counter,5)=y(round(f(j)))-y(j);
  varray(counter,6)=z(round(f(j)))-z(j);
  varray(counter,7:9)=0.;
  counter=counter+1;
end
l_varray=length(varray);
%----------------------------------------------------------------------
%create the periodic array
per_counter=1;
varray2(1:7*l_varray,1:6)=0.;
for i=-1:2:1 ;
  varray2((per_counter-1)*l_varray+1:per_counter*l_varray,1)=varray(:,1)+i*dims(2);
  varray2((per_counter-1)*l_varray+1:per_counter*l_varray,2)=varray(:,2);
  varray2((per_counter-1)*l_varray+1:per_counter*l_varray,3)=varray(:,3);
  varray2((per_counter-1)*l_varray+1:per_counter*l_varray,4)=varray(:,4);
  varray2((per_counter-1)*l_varray+1:per_counter*l_varray,5)=varray(:,5);
  varray2((per_counter-1)*l_varray+1:per_counter*l_varray,6)=varray(:,6);
  per_counter=per_counter+1;
end
for i=-1:2:1 ;
  varray2((per_counter-1)*l_varray+1:per_counter*l_varray,1)=varray(:,1);
  varray2((per_counter-1)*l_varray+1:per_counter*l_varray,2)=varray(:,2)+i*dims(2);
  varray2((per_counter-1)*l_varray+1:per_counter*l_varray,3)=varray(:,3);
  varray2((per_counter-1)*l_varray+1:per_counter*l_varray,4)=varray(:,4);
  varray2((per_counter-1)*l_varray+1:per_counter*l_varray,5)=varray(:,5);
  varray2((per_counter-1)*l_varray+1:per_counter*l_varray,6)=varray(:,6);
  per_counter=per_counter+1;
end ;
for i=-1:2:1 ;
  varray2((per_counter-1)*l_varray+1:per_counter*l_varray,1)=varray(:,1);
  varray2((per_counter-1)*l_varray+1:per_counter*l_varray,2)=varray(:,2);
  varray2((per_counter-1)*l_varray+1:per_counter*l_varray,3)=varray(:,3)+i*dims(2);
  varray2((per_counter-1)*l_varray+1:per_counter*l_varray,4)=varray(:,4);
  varray2((per_counter-1)*l_varray+1:per_counter*l_varray,5)=varray(:,5);
  varray2((per_counter-1)*l_varray+1:per_counter*l_varray,6)=varray(:,6);
  per_counter=per_counter+1;
end
varray2((per_counter-1)*l_varray+1:per_counter*l_varray,1)=varray(:,1);
varray2((per_counter-1)*l_varray+1:per_counter*l_varray,2)=varray(:,2);
varray2((per_counter-1)*l_varray+1:per_counter*l_varray,3)=varray(:,3);
varray2((per_counter-1)*l_varray+1:per_counter*l_varray,4)=varray(:,4);
varray2((per_counter-1)*l_varray+1:per_counter*l_varray,5)=varray(:,5);
varray2((per_counter-1)*l_varray+1:per_counter*l_varray,6)=varray(:,6);
%----------------------------------------------------------------------
% plot3(varray2(:,1),varray2(:,2),varray2(:,3),'r.')
% hold on
% plot3(varray(:,1),varray(:,2),varray(:,3),'b.')
% return
bsize=dims(2);
ns = createns(varray2(:,1:3),'nsmethod','kdtree');
%now we need to define the mesh
%let's get the smoothing length
ts=load('./data/ts.log');
delta=sqrt((bsize^3)/ts(filenumber,6))/2;
for i=1:l_varray ;
  [ind ddist]=knnsearch(ns,varray(i,1:3),'k',200);
  if max(ddist)<2*delta
    disp('i think more points are needed in KD tree')
  end
  for l=1:length(ind)
    factor=SPH_kernel_new(ddist(l),delta);
    %test(counter,1)=ddist(l);
    %test(counter,2)=factor;
    %factor=gaussian_kernel(ddist(l),delta);
    %test(counter,3)=factor;
    %counter=counter+1;
    %factor=SPH_kernel(ddist(l),delta);
    varray(i,7)=varray(i,7)+factor*varray2(ind(l),4);
    varray(i,8)=varray(i,8)+factor*varray2(ind(l),5);
    varray(i,9)=varray(i,9)+factor*varray2(ind(l),6);
  end
  clear ind ddist
end ; 
vort=sqrt(varray(:,7).^2+varray(:,8).^2+varray(:,9).^2);
vort_rms=sqrt(sum(varray(:,7).^2+varray(:,8).^2+varray(:,9).^2)/l_varray);
W=vort;
%return
%now create the povray script %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wmax=max(W); %maximum value for the smoothed vorticity
tones = 1000; %number of different color tones
ColMap=colormap(jet(tones+1));
%tuberad=0.02;
%tuberad=0.02/10;
tuberad=0.0002;
scale = 100.00; %Added in order to get povray working correctly!

tuberad = scale*tuberad;
fid = fopen('mesh.pov','w');
fprintf(fid,'#declare Vortices = \n');
fprintf(fid,'union {\n');
for j=1:l_varray
    dummy_x(1,1)=varray(j,1);
    dummy_x(2,1)=varray(j,4)+varray(j,1);
    dummy_x(1,2)=varray(j,2);
    dummy_x(2,2)=varray(j,5)+varray(j,2);
    dummy_x(1,3)=varray(j,3);
    dummy_x(2,3)=varray(j,6)+varray(j,3);
      wval=W(j); %use when knowing the smoothed circulation
      %if (onlycoherent == 0) | (wval>treshold*wrms) 
        fprintf(fid,'sphere_sweep {\n');
        fprintf(fid,'linear_spline\n');
        N = 2; %linear interpolation
        fprintf(fid,'%5i\n',N);
        dummy_x = scale*dummy_x;
        fprintf(fid,'<%10.4f, %10.4f, %10.4f>, %10.4f\n',[dummy_x(1,1) dummy_x(1,2) dummy_x(1,3) tuberad]);  
        fprintf(fid,'<%10.4f, %10.4f, %10.4f>, %10.4f\n',[dummy_x(2,1) dummy_x(2,2) dummy_x(2,3) tuberad]);
        %fprintf(fid,'tolerance 0.1\n');
        colorInd=round(tones*wval/wmax)+1;    
        rgbvals=ColMap(colorInd,:);
        %fprintf(fid,'pigment { color red 1 green 0 blue 0 }\n');
        fprintf(fid,'pigment { color red %6.4f green %6.4f blue %6.4f }\n',rgbvals);
        fprintf(fid,'finish { ambient 0.2 diffuse 0.99 phong 1 }\n');
        fprintf(fid,'}\n');
      %end %if
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
  
end
function factor=SPH_kernel(r,h)
  q=r/(h);
  if q<=1 
    factor=(2-q)^3-4*(1-q)^3;
  elseif q<=2
    factor=(2-q)^3;
  else
    factor=0.;
  end
  factor=9.97E-4*factor*0.0163/(h^3);
end
function factor=SPH_kernel_new(r,h)
  q=r/(2*h);
  if q<=1.
    factor=1-1.5*(q^2)+0.75*(q^3);
  elseif q<=2
    factor=0.25*((2-q)^3);
  else
    factor=0.;
  end
  factor=9.97E-4*factor/((2*pi*h^2)^(3/2));
end
function factor=gaussian_kernel(r,h)
  factor=exp(-r^2/(2*h^2));
  factor=9.97E-4*factor/((2*pi*h^2)^(3/2));
end
function factor=SPH_kernel_springel(r,h)
  q=r/(2*h);
  if q<=0.5 
    factor=1-6*(q^2)+6*(q^3);
  elseif q<=1
    factor=2*((1-q)^3);
  else
    factor=0.;
  end
  factor=9.97E-4*2.54647909*factor/(h^3);
end