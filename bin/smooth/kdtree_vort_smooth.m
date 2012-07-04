function [L_inside L_outside L_total vort_rms]=kdtree_vort_smooth(filenumber,plot,eval_rms)
global dims
global x y z
global f
global number_of_particles
global vort total_length vort_rms
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
total_length=sqrt(varray(:,4).^2+varray(:,5).^2+varray(:,6).^2);
%return
L_total=sum(total_length);
for i=1:10
  eval_rms2=eval_rms+0.12*(i-1)
  ind=vort>eval_rms2*vort_rms;
  L_inside(i)=sum(total_length(ind));
  L_outside(i)=L_total-L_inside(i);
  clear ind
end
if plot==1
  vort_norm=vort-min(vort);
  rainbow_scale=299/max(vort_norm);
  vort_norm=vort_norm*rainbow_scale;
  store_caxis=([min(vort) max(vort)]);
  rainbowcmap=colormap(hot(400));
  for i=1:l_varray
%    plot3([varray(i,1) (varray(i,4)+varray(i,1))],[varray(i,2) (varray(i,5)+varray(i,2))],[varray(i,3) (varray(i,6)+varray(i,3))],'-','Color',rainbowcmap(max(1,ceil(vort_norm(i))),:),'LineWidth',1.)
    subplot(2,2,2)
    if vort(i)>1*vort_rms
          plot3([varray(i,1) (varray(i,4)+varray(i,1))],...
          [varray(i,2) (varray(i,5)+varray(i,2))],...
          [varray(i,3) (varray(i,6)+varray(i,3))],...
          '-','Color',rainbowcmap(max(1,ceil(10*vort(i))),:),...
          'LineWidth',1.)

    hold on
    end
    subplot(2,2,4)
    if vort(i)<1*vort_rms
          plot3([varray(i,1) (varray(i,4)+varray(i,1))],...
          [varray(i,2) (varray(i,5)+varray(i,2))],...
          [varray(i,3) (varray(i,6)+varray(i,3))],...
          '-','Color',rainbowcmap(max(1,ceil(10*vort(i))),:),...
          'LineWidth',1.)

    hold on
    end
    subplot(2,2,[1 3])
    plot3([varray(i,1) (varray(i,4)+varray(i,1))],...
          [varray(i,2) (varray(i,5)+varray(i,2))],...
          [varray(i,3) (varray(i,6)+varray(i,3))],...
          '-','Color',rainbowcmap(max(1,ceil(10*vort(i))),:),...
          'LineWidth',1.)

    hold on
  end
  subplot(2,2,2)
  camproj('perspective')
  box on
  axis([-dims(2)/2 dims(2)/2 -dims(2)/(2*dims(7)) dims(2)/(2*dims(7)) -dims(2)/(2*dims(7)) dims(2)/(2*dims(7))]);
  daspect([1 dims(7) dims(7)])
  set(gca,'xtick',[])
  set(gca,'ytick',[])
  set(gca,'ztick',[])
  subplot(2,2,4)
  camproj('perspective')
  box on
  axis([-dims(2)/2 dims(2)/2 -dims(2)/(2*dims(7)) dims(2)/(2*dims(7)) -dims(2)/(2*dims(7)) dims(2)/(2*dims(7))]);
  daspect([1 dims(7) dims(7)])
  set(gca,'xtick',[])
  set(gca,'ytick',[])
  set(gca,'ztick',[])
  subplot(2,2,[1 3])
  camproj('perspective')
  box on
  axis([-dims(2)/2 dims(2)/2 -dims(2)/(2*dims(7)) dims(2)/(2*dims(7)) -dims(2)/(2*dims(7)) dims(2)/(2*dims(7))]);
  daspect([1 dims(7) dims(7)])
  set(gca,'xtick',[])
  set(gca,'ytick',[])
  set(gca,'ztick',[])
  caxis([0 1])
  colorbar
  if 1==1
    whitebg('k')
    set(gcf,'InvertHardcopy','off');
  end
  hold off
end
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