function [L_inside L_outside L_total rms_w]=kdtree_smooth(filenumber,msize,plot,eval_rms)
global dims
global x y z
global f
global number_of_particles
if nargin<3
  plot=1;
  eval_rms=2.;
elseif nargin<4
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
  counter=counter+1;
end
l_varray=length(varray);
counter=2;
bsize=dims(2);
% for i=-1:1 ; for j=-1:1 ; for k=-1:1
%   if i==0 && j==0 && k==0 ;  continue ; end
%   varray((counter-1)*l_varray+1:counter*l_varray,1)=varray(1:l_varray,1)+i*bsize;
%   varray((counter-1)*l_varray+1:counter*l_varray,2)=varray(1:l_varray,2)+j*bsize;
%   varray((counter-1)*l_varray+1:counter*l_varray,3)=varray(1:l_varray,3)+k*bsize;
%   counter=counter+1;
% end ;  end ; end ;
%plot3(varray(:,1),varray(:,2),varray(:,3),'.')
%return
ns = createns(varray(:,1:3),'nsmethod','kdtree');
%now we need to define the mesh
res=(bsize/msize);
clear wmesh
wmesh(1:msize,1:msize,1:msize,1:3)=0.;
%let's get the smoothing length
ts=load('./data/ts.log');
delta=sqrt((bsize^3)/ts(filenumber,6))/2;
counter=1;
for i=1:msize ;
  i ;
  for j=1:msize ; for k=1:msize ;
  xx(1)= res*(2*i-1)/2.-(bsize/2.);
  xx(2)= res*(2*j-1)/2.-(bsize/2.);
  xx(3)= res*(2*k-1)/2.-(bsize/2.);
  [ind ddist]=knnsearch(ns,xx,'k',400);
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
    wmesh(k,j,i,1)=wmesh(k,j,i,1)+factor*varray(ind(l),4);
    wmesh(k,j,i,2)=wmesh(k,j,i,2)+factor*varray(ind(l),5);
    wmesh(k,j,i,3)=wmesh(k,j,i,3)+factor*varray(ind(l),6);
  end
  clear ind ddist
end ; end ; end ; 
if plot==1
  mesh_iso(linspace(-bsize/2,bsize/2,msize),wmesh(:,:,:,1),wmesh(:,:,:,2),wmesh(:,:,:,3),msize,'smooth')
end
%---------------now work out length inside and outside structures----
mod_w=(wmesh(:,:,:,1).^2+wmesh(:,:,:,2).^2+wmesh(:,:,:,3).^2);
rms_w=sqrt(sum(sum(sum(mod_w)))/(msize^3));
[sw3,sw2,sw1]=ind2sub(size(mod_w),find(sqrt(mod_w)>eval_rms*rms_w));
for i=1:msize
  xmesh(i)=((2*i-1)/2)*bsize/msize-bsize/2;
end
L_inside=0.;
L_total=0.;
parfor l=1:length(varray)
  L_total=L_total+sqrt(varray(l,4)^2+varray(l,5)^2+varray(l,6)^2);
  hold on
  for i=1:length(sw1)
      if sw1(i)==msize || sw2(i)==msize || sw3(i)==msize
        continue
      end
      if varray(l,1)>xmesh(sw1(i)) && varray(l,1)<xmesh(sw1(i)+1) &&  varray(l,2)>xmesh(sw2(i)) && varray(l,2)<xmesh(sw2(i)+1) && varray(l,3)>xmesh(sw3(i)) && varray(l,3)<xmesh(sw3(i)+1)
        L_inside=L_inside+sqrt(varray(l,4)^2+varray(l,5)^2+varray(l,6)^2);
        break
      end        
    end
end
L_outside=L_total-L_inside;
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