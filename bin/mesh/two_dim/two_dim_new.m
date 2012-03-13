function two_dim_slice(filenumber,fit,fitnorm,fitextra,firecolor)
close all
if nargin<2
  do_fit=0;
  firecolor=0;
else
  do_fit=1;
end
if nargin<3
  if do_fit==1
    fitnorm=fit;
  end
  firecolor=0;
end
if nargin<4
  fitextra=[];
  firecolor=0;
end
if nargin<5
  firecolor=0;
end
filename=sprintf('data/vel_slice_2D%04d.dat',filenumber);
fid=fopen(filename);
if fid<0
  disp('2D slice file does not exist, exiting script')
  return
end
if firecolor==1
  disp('using fireprint colour scheme')
end
load data/dims.log;
msize=dims(8);
if (msize==0) 
  disp('2D mesh size is zero exiting script')
  return
end
%t=fread(fid,1,'float64');
x=fread(fid,msize,'float64');
unormx=fread(fid,msize^2,'float64');
unormy=fread(fid,msize^2,'float64');
unormz=fread(fid,msize^2,'float64');

usupx=fread(fid,msize^2,'float64');
usupy=fread(fid,msize^2,'float64');
usupz=fread(fid,msize^2,'float64');

unormx=reshape(unormx,msize,msize);
unormy=reshape(unormy,msize,msize);
unormz=reshape(unormz,msize,msize);
usupx=reshape(usupx,msize,msize);
usupy=reshape(usupy,msize,msize);
usupz=reshape(usupz,msize,msize);
xx=x;
yy=x;
figure('Name','2D (super) velocity slice')
subplot(2,2,1)
imagesc(xx,yy,usupx) ; shading interp
set(gca,'YDir','normal')
axis square
set(gca,'FontSize',14)
xlabel('u_x','FontSize',14)
colorbar
subplot(2,2,2)
imagesc(xx,yy,usupy) ; shading interp
set(gca,'YDir','normal')
axis square
set(gca,'FontSize',14)
xlabel('u_y','FontSize',14)
colorbar
subplot(2,2,3)
imagesc(xx,yy,usupz) ; shading interp
set(gca,'YDir','normal')
axis square
set(gca,'FontSize',14)
xlabel('u_z','FontSize',14)
colorbar
subplot(2,2,4)
imagesc(xx,yy,log(sqrt(usupx.^2+usupy.^2+usupz.^2))) ; shading interp
set(gca,'YDir','normal')
axis square
set(gca,'FontSize',14)
xlabel('log|u|','FontSize',14)
if firecolor==1
  colormap(fireprint)
end
colorbar

figure('Name','2D (super) |u| slice')
uu=sqrt(usupx.^2+usupy.^2+usupz.^2);
vcoff=1. ;
index = find(uu > vcoff);
uu(index) = vcoff;
clear index
imagesc(xx,yy,uu) ; shading interp
 set(gca,'YDir','normal')
set(gca,'FontSize',14)
xlabel('|u|','FontSize',14)
if firecolor==1
  colormap(fireprint)
end
colorbar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if std(std(sqrt(unormx.^2+unormy.^2+unormz.^2)))>0
  figure('Name','2D (normal) velocity slice')
  subplot(2,2,1)
  imagesc(xx,yy,unormx) ; shading interp
  set(gca,'YDir','normal')
  axis square
  set(gca,'FontSize',14)
  xlabel('u_x','FontSize',14)
  colorbar
  subplot(2,2,2)
  imagesc(xx,yy,unormy) ; shading interp
  set(gca,'YDir','normal')
  axis square
  set(gca,'FontSize',14)
  xlabel('u_y','FontSize',14)
  colorbar
  subplot(2,2,3)
  imagesc(xx,yy,unormz) ; shading interp
  set(gca,'YDir','normal')
  axis square
  set(gca,'FontSize',14)
  xlabel('u_z','FontSize',14)
  colorbar
  subplot(2,2,4)
  imagesc(xx,yy,sqrt(unormx.^2+unormy.^2+unormz.^2)) ; shading interp
  set(gca,'YDir','normal')
  axis square
  set(gca,'FontSize',14)
  xlabel('|u|','FontSize',14)
  if firecolor==1
    colormap(fireprint)
  end
  colorbar

  figure('Name','2D (normal) |u| slice')
  imagesc(xx,yy,sqrt(unormx.^2+unormy.^2+unormz.^2)) ; shading interp
  set(gca,'YDir','normal')
  set(gca,'FontSize',14)
  xlabel('x','FontSize',14)
  ylabel('y','FontSize',14)

  if firecolor==1
    colormap(fireprint)
  end
  colorbar
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n=msize;
fux=fftn(usupx)/(n^2);
fuy=fftn(usupy)/(n^2);
fuz=fftn(usupz)/(n^2);
energyr=real(fux).^2+real(fuy).^2+real(fuz).^2;
energyi=imag(fux).^2+imag(fuy).^2+imag(fuz).^2;
midpt=n/2+1;
spect(1:1.5*n)=0.;
for i=1:n
  for j=1:n
    ii=i;
    jj=j;
    if ii>midpt ; ii=n-ii+1;  end ;
    if jj>midpt ; jj=n-jj+1;  end ;
    r=int16(sqrt(ii^2+jj^2));
    spect(r)=spect(r)+energyr(i,j)+energyi(i,j);
  end
end
cutoff=100;
figure('Name','Superfluid energy spectrum')
k=1:midpt;
k2=(floor(midpt/2):midpt);
loglog(k(1:midpt-cutoff),spect(1:midpt-cutoff),'LineWidth',2)
xlabel('log k','FontSize',14) ; ylabel('log E(k)','FontSize',14)
axis tight
set(gca,'FontSize',14)
if do_fit==1
  disp(sprintf('fitting a slope of %f to plot',fit))
  dummy_spect=k.^(fit);
  scaling_factor=sum(spect(1:10))/sum(dummy_spect(1:10));
  dummy_spect=dummy_spect*scaling_factor;
  hold on
  loglog(k,dummy_spect,'--k','LineWidth',2)
  if isempty(fitextra)==0
    disp(sprintf('fitting a slope of %f to plot',fitextra))
    dummy_spect2=k2.^(fitextra);
    scaling_factor=sum(spect(floor(midpt/2)+20))/sum(dummy_spect2(20));
    %scaling_factor=sum(spect(1:length(spect)))/sum(dummy_spect((1:length(dummy_spect))));
    dummy_spect2=dummy_spect2*scaling_factor;
    loglog(k2,dummy_spect2,'--r','LineWidth',2)
  end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Name','Superfluid (-5/3) compensated energy spectrum')
comp_spec=(k.^(5/3));
loglog(k(1:midpt-cutoff),comp_spec(1:midpt-cutoff).*spect(1:midpt-cutoff),'LineWidth',2)
xlabel('log k','FontSize',14) ; ylabel('log E(k)*k^{5/3}','FontSize',14)
axis tight
set(gca,'FontSize',14)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Name','Superfluid (-1) compensated energy spectrum')
comp_spec=(k.^1);
loglog(k(1:midpt-cutoff),comp_spec(1:midpt-cutoff).*spect(1:midpt-cutoff),'LineWidth',2)
xlabel('log k','FontSize',14) ; ylabel('log E(k)*k','FontSize',14)
axis tight
set(gca,'FontSize',14)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if std(std(sqrt(unormx.^2+unormy.^2+unormz.^2)))>0
  fux=fftn(unormx)/(n^2);
  fuy=fftn(unormy)/(n^2);
  fuz=fftn(unormz)/(n^2);
  energyr=real(fux).^2+real(fuy).^2+real(fuz).^2;
  energyi=imag(fux).^2+imag(fuy).^2+imag(fuz).^2;
  spectnorm(1:1.5*n)=0.;
  for i=1:n
    for j=1:n
      ii=i;
      jj=j;
      if ii>midpt ; ii=n-ii+1; end ;
      if jj>midpt ; jj=n-jj+1; end ;
      r=int16(sqrt(ii^2+jj^2));
      spectnorm(r)=spectnorm(r)+energyr(i,j)+energyi(i,j);
    end
  end
  figure('Name','Normal fluid energy spectrum')
  loglog(k,spectnorm(1:midpt),'LineWidth',2)
  xlabel('log k','FontSize',14) ; ylabel('log E(k)','FontSize',14)
  axis tight
  set(gca,'FontSize',14)
  if do_fit==1
    disp(sprintf('fitting a slope of %f to plot',fitnorm))
    dummy_spect=k.^(fitnorm);
    scaling_factor=sum(spectnorm(1:10))/sum(dummy_spect(1:10));
    dummy_spect=dummy_spect*scaling_factor;
    hold on
    loglog(k,dummy_spect,'--k','LineWidth',2)
  end
end
