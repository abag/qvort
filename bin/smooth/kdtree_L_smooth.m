function kdtree_vort_smooth(filenumber,msize)
global dims
global x y z
global f
global number_of_particles
global vort total_length vort_rms
dims=load('./data/dims.log');
%global test
filename=sprintf('data/var%04d.log',filenumber);
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
counter=1;
for j=1:number_of_particles
  if f(j)==0 ; continue ; end
  dist=sqrt((x(j)-x(round(f(j))))^2+(y(j)-y(round(f(j))))^2+(z(j)-z(round(f(j))))^2);
  if (dist>0.5*dims(2)) ; continue ; end
  varray(counter,1)=x(j);
  varray(counter,2)=y(j);
  varray(counter,3)=z(j);
  varray(counter,4)=sqrt((x(round(f(j)))-x(j))^2+(x(round(f(j)))-x(j))^2+(x(round(f(j)))-x(j))^2);
  counter=counter+1;
end
l_varray=length(varray);
%----------------------------------------------------------------------
%create the periodic array
per_counter=1;
varray2(1:5*l_varray,1:4)=0.;
for i=-1:2:1 ;
  varray2((per_counter-1)*l_varray+1:per_counter*l_varray,1)=varray(:,1)+i*dims(2);
  varray2((per_counter-1)*l_varray+1:per_counter*l_varray,2)=varray(:,2);
  varray2((per_counter-1)*l_varray+1:per_counter*l_varray,3)=varray(:,3);
  varray2((per_counter-1)*l_varray+1:per_counter*l_varray,4)=varray(:,4);
  per_counter=per_counter+1;
end
for i=-1:2:1 ;
  varray2((per_counter-1)*l_varray+1:per_counter*l_varray,1)=varray(:,1);
  varray2((per_counter-1)*l_varray+1:per_counter*l_varray,2)=varray(:,2)+i*dims(2);
  varray2((per_counter-1)*l_varray+1:per_counter*l_varray,3)=varray(:,3);
  varray2((per_counter-1)*l_varray+1:per_counter*l_varray,4)=varray(:,4);
  per_counter=per_counter+1;
end ;
for i=-1:2:1 ;
  varray2((per_counter-1)*l_varray+1:per_counter*l_varray,1)=varray(:,1);
  varray2((per_counter-1)*l_varray+1:per_counter*l_varray,2)=varray(:,2);
  varray2((per_counter-1)*l_varray+1:per_counter*l_varray,3)=varray(:,3)+i*dims(2);
  varray2((per_counter-1)*l_varray+1:per_counter*l_varray,4)=varray(:,4);
  per_counter=per_counter+1;
end
varray2((per_counter-1)*l_varray+1:per_counter*l_varray,1)=varray(:,1);
varray2((per_counter-1)*l_varray+1:per_counter*l_varray,2)=varray(:,2);
varray2((per_counter-1)*l_varray+1:per_counter*l_varray,3)=varray(:,3);
varray2((per_counter-1)*l_varray+1:per_counter*l_varray,4)=varray(:,4);
%----------------------------------------------------------------------
 %plot3(varray2(:,1),varray2(:,2),varray2(:,3),'r.')
 %hold on
 %plot3(varray(:,1),varray(:,2),varray(:,3),'b.')
 %return
 bsize=dims(2);
 res=(bsize/msize);
 clear Lmesh
 Lmesh(1:msize,1:msize,1:msize)=0.;
 ns = createns(varray2(:,1:3),'nsmethod','kdtree');
 %now we need to define the mesh
 %let's get the smoothing length
 ts=load('./data/ts.log');
 delta=sqrt((bsize^3)/ts(filenumber,6));
for i=1:msize ;
  i 
  for j=1:msize ; for k=1:msize ;
  xx(1)= res*(2*i-1)/2.-(bsize/2.);
  xx(2)= res*(2*j-1)/2.-(bsize/2.);
  xx(3)= res*(2*k-1)/2.-(bsize/2.);
  [ind ddist]=knnsearch(ns,xx,'k',1400);
  if max(ddist)<2*delta
    disp('i think more points are needed in KD tree')
  end
  for l=1:length(ind)
    factor=SPH_kernel_new(ddist(l),delta);
    Lmesh(k,j,i)=Lmesh(k,j,i,1)+factor*varray2(ind(l),4);
  end
  clear ind ddist
end ; end ; end ; 
xmesh=linspace(-bsize/2,bsize/2,msize);
if 1==1
  isosurface(Lmesh)
end
save out.mat Lmesh xmesh
 function factor=SPH_kernel_new(r,h)
   %q=r/(2*h);
   q=r/h;
   if q<=1.
     factor=1-1.5*(q^2)+0.75*(q^3);
   elseif q<=2
     factor=0.25*((2-q)^3);
   else
     factor=0.;
   end
