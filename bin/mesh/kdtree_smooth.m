function kdtree_smooth(filenumber,msize)
global dims
global x y z
global f
global number_of_particles
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
  counter=counter+1;
end
l_varray=length(varray);
counter=2;
bsize=dims(2);
for i=-1:1 ; for j=-1:1 ; for k=-1:1
  if i==0 && j==0 && k==0 ;  continue ; end
  varray((counter-1)*l_varray+1:counter*l_varray,1)=varray(1:l_varray,1)+i*bsize;
  varray((counter-1)*l_varray+1:counter*l_varray,2)=varray(1:l_varray,2)+j*bsize;
  varray((counter-1)*l_varray+1:counter*l_varray,3)=varray(1:l_varray,3)+k*bsize;
  counter=counter+1;
end ;  end ; end ;
%plot3(varray(:,1),varray(:,2),varray(:,3),'.')
%return
ns = createns(varray(:,1:3),'nsmethod','kdtree');
%now we need to define the mesh
res=(bsize/msize);
wmesh(1:msize,1:msize,1:msize,1:3)=0.;
%let's get the smoothing length
ts=load('./data/ts.log');
delta=0.01 ; %1/sqrt(2*ts(filenumber,6)/(bsize^3));
for i=1:msize ;
  i 
  for j=1:msize ; for k=1:msize ;
  xx(1)= res*(2*i-1)/2.-(bsize/2.);
  xx(2)= res*(2*j-1)/2.-(bsize/2.);
  xx(3)= res*(2*k-1)/2.-(bsize/2.);
  [ind ddist]=knnsearch(ns,xx,'k',200);
  for l=1:length(ind)
    factor=SPH_kernel(ddist(l),delta);
    wmesh(k,j,i,1)=wmesh(k,j,i,1)+factor*varray(ind(l),4);
    wmesh(k,j,i,2)=wmesh(k,j,i,2)+factor*varray(ind(l),5);
    wmesh(k,j,i,3)=wmesh(k,j,i,3)+factor*varray(ind(l),6);
  end
  clear ind ddist
end ; end ; end ; 
mesh_iso(linspace(-bsize/2,bsize/2,msize),wmesh(:,:,:,1),wmesh(:,:,:,2),wmesh(:,:,:,3),msize,'smooth')
end
function factor=SPH_kernel(r,h)
  q=r/(2*h);
  if q<=0.5 
    factor=1-6*(q^2)+6*(q^3);
  elseif q<=1
    factor=2*((1-q)^3);
  else
    factor=0.;
  end
end