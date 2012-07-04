function [L_inside L_outside L_total vort_rms]=kdtree_spec(filenumber,plot)
global dims
global x y z
global f
global number_of_particles
global meshu meshu2 meshu3 
global meshu4 meshu5 meshu6 
global meshu7 meshu8 meshu9
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
  if (dist>0.5*dims(2))
    %periodic BC's
    varray(counter,4)=x(round(f(j)))-x(j);
    if varray(counter,4)>0.5*dims(2)
      varray(counter,4)=dims(2)-varray(counter,4);
    end
    varray(counter,5)=y(round(f(j)))-y(j);
    if varray(counter,5)>0.5*dims(2)
      varray(counter,5)=dims(2)-varray(counter,5);
    end
    varray(counter,6)=z(round(f(j)))-z(j);
    if varray(counter,6)>0.5*dims(2)
      varray(counter,6)=dims(2)-varray(counter,6);
    end
  else
    varray(counter,4)=x(round(f(j)))-x(j);
    varray(counter,5)=y(round(f(j)))-y(j);
    varray(counter,6)=z(round(f(j)))-z(j);
  end
  varray(counter,7:9)=0.;
  counter=counter+1;
end
l_varray=length(varray);
%----------------------------------------------------------------------
%create the periodic array
per_counter=1;
varray2(1:27*l_varray,1:9)=0.;
for i=-1:1 ; for j=-1:1 ; for k=-1:1
  varray2((per_counter-1)*l_varray+1:per_counter*l_varray,1)=varray(:,1)+i*dims(2);
  varray2((per_counter-1)*l_varray+1:per_counter*l_varray,2)=varray(:,2)+j*dims(2);
  varray2((per_counter-1)*l_varray+1:per_counter*l_varray,3)=varray(:,3)+k*dims(2);
  varray2((per_counter-1)*l_varray+1:per_counter*l_varray,4)=varray(:,4);
  varray2((per_counter-1)*l_varray+1:per_counter*l_varray,5)=varray(:,5);
  varray2((per_counter-1)*l_varray+1:per_counter*l_varray,6)=varray(:,6);
  per_counter=per_counter+1;
end ; end ; end
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
%----------now make vort 
vort_rms=sqrt(sum(varray(:,7).^2+varray(:,8).^2+varray(:,9).^2)/l_varray);
%----------now make vort 26 times longer for periodicity--------------
%create the periodic array
per_counter=1;
vort2(1:27*l_varray)=0.;
for i=-1:1 ; for j=-1:1 ; for k=-1:1
  vort2((per_counter-1)*l_varray+1:per_counter*l_varray)=vort(:);
  per_counter=per_counter+1;
end ; end ; end
%---------------------------------------------------------------------

if plot==1
  vort_norm=vort-min(vort);
  rainbow_scale=299/max(vort_norm);
  vort_norm=vort_norm*rainbow_scale;
  store_caxis=([min(vort) max(vort)]);
  rainbowcmap=colormap(hot(400));
  for i=1:l_varray
    plot3([varray(i,1) (varray(i,4)+varray(i,1))],...
          [varray(i,2) (varray(i,5)+varray(i,2))],...
          [varray(i,3) (varray(i,6)+varray(i,3))],...
          '-','Color',rainbowcmap(max(1,ceil(vort_norm(i))),:),'LineWidth',1.)
    hold on
  end
  camproj('perspective')
  box on
  axis([-dims(2)/2 dims(2)/2 -dims(2)/(2*dims(7)) dims(2)/(2*dims(7)) -dims(2)/(2*dims(7)) dims(2)/(2*dims(7))]);
  daspect([1 dims(7) dims(7)])
  set(gca,'xtick',[])
  set(gca,'ytick',[])
  set(gca,'ztick',[])
  caxis(store_caxis)
  colorbar
  if 1==1
    whitebg('k')
    set(gcf,'InvertHardcopy','off');
  end
  hold off
end
%-----------------------------now compute spectra--------------------
nmesh=256;
npoints=500;
for i=1:nmesh
  xmesh(i)=-dims(2)/2+((2*i-1)/(2*nmesh))*dims(2);
end
%-------------------------------------------------------------------
meshu(1:nmesh,1:nmesh,1:3)=0. ;
meshu2(1:nmesh,1:nmesh,1:3)=0. ;
meshu3(1:nmesh,1:nmesh,1:3)=0. ;
meshu4(1:nmesh,1:nmesh,1:3)=0. ;
meshu5(1:nmesh,1:nmesh,1:3)=0. ;
meshu6(1:nmesh,1:nmesh,1:3)=0. ;
meshu7(1:nmesh,1:nmesh,1:3)=0. ;
meshu8(1:nmesh,1:nmesh,1:3)=0. ;
meshu9(1:nmesh,1:nmesh,1:3)=0. ;
for i=1:nmesh
  i
  for j=1:nmesh
    %get the location and distance to the nearest points
    dummy_loc=[xmesh(i) xmesh(j) 0.];
    [ind ddist]=knnsearch(ns,dummy_loc,'k',npoints);
    %loop over these points
    for k=1:npoints
      a_bs=ddist(k)^2;
      b_bs=2*dot((varray2(ind(k),1:3)-dummy_loc),varray2(ind(k),4:6));
      c_bs=norm(varray2(ind(k),4:6))^2;
      if ((4*a_bs*c_bs-b_bs^22)==0) ; continue ; end  ;%avoid 1/0!
      u_bs=cross((varray2(ind(k),1:3)-dummy_loc),varray2(ind(k),4:6));
      u_bs=u_bs*9.97E-4/((2*pi)*(4*a_bs*c_bs-b_bs^2));
      u_bs=u_bs*((2*c_bs+b_bs)/sqrt(a_bs+b_bs+c_bs)-(b_bs/sqrt(a_bs)));
      meshu(j,i,1:3)=squeeze(meshu(j,i,1:3))+u_bs(1:3)';
      if vort2(ind(k))>1.2*vort_rms
          meshu9(j,i,1:3)=squeeze(meshu9(j,i,1:3))+u_bs(1:3)';
        if vort2(ind(k))>1.4*vort_rms
          meshu2(j,i,1:3)=squeeze(meshu2(j,i,1:3))+u_bs(1:3)';
          if vort2(ind(k))>1.7*vort_rms
            meshu3(j,i,1:3)=squeeze(meshu3(j,i,1:3))+u_bs(1:3)';
            if vort2(ind(k))>2*vort_rms
              meshu4(j,i,1:3)=squeeze(meshu4(j,i,1:3))+u_bs(1:3)';
            end
          end
        end
      end
      if vort2(ind(k))<1.7*vort_rms
        meshu5(j,i,1:3)=squeeze(meshu5(j,i,1:3))+u_bs(1:3)';
        if vort2(ind(k))<1.4*vort_rms
          meshu6(j,i,1:3)=squeeze(meshu6(j,i,1:3))+u_bs(1:3)';
          if vort2(ind(k))<1.2*vort_rms
            meshu7(j,i,1:3)=squeeze(meshu7(j,i,1:3))+u_bs(1:3)';
            if vort2(ind(k))<1*vort_rms
              meshu8(j,i,1:3)=squeeze(meshu8(j,i,1:3))+u_bs(1:3)';
            end
          end
        end
      end
    end
  end
end
pcolor(sqrt(meshu(:,:,1).^2+meshu(:,:,2).^2+meshu(:,:,3).^2)) ; shading interp
%------------------spectra------------------
fux=fftn(meshu(:,:,1))/(nmesh^2);
fuy=fftn(meshu(:,:,2))/(nmesh^2);
fuz=fftn(meshu(:,:,3))/(nmesh^2);
energyr=real(fux).^2+real(fuy).^2+real(fuz).^2;
energyi=imag(fux).^2+imag(fuy).^2+imag(fuz).^2;
%-------------------------------------
fux2=fftn(meshu2(:,:,1))/(nmesh^2);
fuy2=fftn(meshu2(:,:,2))/(nmesh^2);
fuz2=fftn(meshu2(:,:,3))/(nmesh^2);
energyr2=real(fux2).^2+real(fuy2).^2+real(fuz2).^2;
energyi2=imag(fux2).^2+imag(fuy2).^2+imag(fuz2).^2;
%-------------------------------------
fux3=fftn(meshu3(:,:,1))/(nmesh^2);
fuy3=fftn(meshu3(:,:,2))/(nmesh^2);
fuz3=fftn(meshu3(:,:,3))/(nmesh^2);
energyr3=real(fux3).^2+real(fuy3).^2+real(fuz3).^2;
energyi3=imag(fux3).^2+imag(fuy3).^2+imag(fuz3).^2;
%-------------------------------------
fux5=fftn(meshu5(:,:,1))/(nmesh^2);
fuy5=fftn(meshu5(:,:,2))/(nmesh^2);
fuz5=fftn(meshu5(:,:,3))/(nmesh^2);
energyr5=real(fux5).^2+real(fuy5).^2+real(fuz5).^2;
energyi5=imag(fux5).^2+imag(fuy5).^2+imag(fuz5).^2;
%-------------------------------------
fux6=fftn(meshu6(:,:,1))/(nmesh^2);
fuy6=fftn(meshu6(:,:,2))/(nmesh^2);
fuz6=fftn(meshu6(:,:,3))/(nmesh^2);
energyr6=real(fux6).^2+real(fuy6).^2+real(fuz6).^2;
energyi6=imag(fux6).^2+imag(fuy6).^2+imag(fuz6).^2;
%-------------------------------------
fux7=fftn(meshu7(:,:,1))/(nmesh^2);
fuy7=fftn(meshu7(:,:,2))/(nmesh^2);
fuz7=fftn(meshu7(:,:,3))/(nmesh^2);
energyr7=real(fux7).^2+real(fuy7).^2+real(fuz7).^2;
energyi7=imag(fux7).^2+imag(fuy7).^2+imag(fuz7).^2;
%-------------------------------------
fux8=fftn(meshu8(:,:,1))/(nmesh^2);
fuy8=fftn(meshu8(:,:,2))/(nmesh^2);
fuz8=fftn(meshu8(:,:,3))/(nmesh^2);
energyr8=real(fux8).^2+real(fuy8).^2+real(fuz8).^2;
energyi8=imag(fux8).^2+imag(fuy8).^2+imag(fuz8).^2;
%-------------------------------------
midpt=nmesh/2+1;
spect(1:1.5*nmesh)=0.;
spect2(1:1.5*nmesh)=0.;
spect3(1:1.5*nmesh)=0.;
spect5(1:1.5*nmesh)=0.;
spect6(1:1.5*nmesh)=0.;
spect7(1:1.5*nmesh)=0.;
spect8(1:1.5*nmesh)=0.;
whitebg('w')
figure('Name','E(k)')
for i=1:nmesh
  for j=1:nmesh
    ii=i;
    jj=j;
    if ii>midpt ; ii=nmesh-ii+1; end ;
    if jj>midpt ; jj=nmesh-jj+1; end ;
    r=int16(sqrt(ii^2+jj^2));
    spect(r)=spect(r)+energyr(i,j)+energyi(i,j);
    spect2(r)=spect2(r)+energyr2(i,j)+energyi2(i,j);
    spect3(r)=spect3(r)+energyr3(i,j)+energyi3(i,j);
    spect5(r)=spect5(r)+energyr5(i,j)+energyi5(i,j);
    spect6(r)=spect6(r)+energyr6(i,j)+energyi6(i,j);
    spect7(r)=spect7(r)+energyr7(i,j)+energyi7(i,j);
    spect8(r)=spect8(r)+energyr8(i,j)+energyi8(i,j);
  end
end
k=(1:midpt)*(2*pi/dims(2));
loglog(k(1:midpt-10),spect(1:midpt-10),'Color','k','LineWidth',1.5)
hold on
loglog(k(1:midpt-10),spect2(1:midpt-10),'Color','r','LineWidth',1.5)
loglog(k(1:midpt-10),spect3(1:midpt-10),'Color','b','LineWidth',1.5)
loglog(k(1:midpt-10),spect5(1:midpt-10),'Color','g','LineWidth',1.5)
loglog(k(1:midpt-10),spect6(1:midpt-10),'Color','m','LineWidth',1.5)
loglog(k(1:midpt-10),spect7(1:midpt-10),'Color','m','LineWidth',1.5)
loglog(k(1:midpt-10),spect8(1:midpt-10),'Color','m','LineWidth',1.5)
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