function two_dim_octave(filenumber,fit)
close all
if nargin<2
  do_fit=0;
else
  do_fit=1;
end
%%%%%%%%%%LOAD IN DATA%%%%%%%%%
filename=sprintf('./data/vel_slice_2D%04d.dat',filenumber);
fid=fopen(filename);
if fid<0
  disp('2D slice file does not exist, exiting script')
  return
end
A=fread(fid,'float64');
s=length(A); s=s/8; s=sqrt(s);
B=reshape(A,8,s,s);
x=squeeze(B(1,:,:));
y=squeeze(B(2,:,:));
xx=squeeze(x(1,:));
yy=squeeze(y(:,1));
usupx=squeeze(B(3,:,:));
usupy=squeeze(B(4,:,:));
usupz=squeeze(B(5,:,:));
unormx=squeeze(B(6,:,:));
unormy=squeeze(B(6,:,:));
unormz=squeeze(B(6,:,:));
%%%%%%%%%%%%%%%%%%SUPER SLICE%%%%%%%%%%%%%
figure('Name','2D (super) velocity slice')
subplot(2,2,1)
imagesc(xx(1,:),yy(:,1),usupx(:,:)) ; shading interp
axis square
set(gca,'FontSize',14)
xlabel('u_x','FontSize',14)
colorbar
subplot(2,2,2)
imagesc(xx(1,:),yy(:,1),usupy(:,:)) ; shading interp
axis square
set(gca,'FontSize',14)
xlabel('u_y','FontSize',14)
colorbar
subplot(2,2,3)
imagesc(xx(1,:),yy(:,1),usupz(:,:)) ; shading interp
axis square
set(gca,'FontSize',14)
xlabel('u_z','FontSize',14)
colorbar
subplot(2,2,4)
uu=sqrt(usupx.^2+usupy.^2+usupz.^2);
imagesc(xx(1,:),yy(:,1),uu(:,:)) ; shading interp
axis square
set(gca,'FontSize',14)
xlabel('log|u|','FontSize',14)
colorbar
%%%%%%%%%%%%%%%%%%SUPER SLICE2%%%%%%%%%%%%%
figure('Name','2D (super) log|u| slice')
imagesc(xx(1,:),yy(:,1),log(uu(:,:))) ; shading interp
set(gca,'FontSize',14)
xlabel('|u|','FontSize',14)
colorbar
%%%%%%%%%%%%%%%%%NORMALSLICE%%%%%%%%%%%%%%%
if std(sqrt(unormx.^2+unormy.^2+unormz.^2))>0
  figure('Name','2D (normal) velocity slice')
  subplot(2,2,1)
  imagesc(xx(1,:),yy(:,1),unormx(:,:)) ; shading interp
  axis square
  set(gca,'FontSize',14)
  xlabel('u_x','FontSize',14)
  colorbar
  subplot(2,2,2)
  imagesc(xx(1,:),yy(:,1),unormy(:,:)) ; shading interp
  axis square
  set(gca,'FontSize',14)
  xlabel('u_y','FontSize',14)
  colorbar
  subplot(2,2,3)
  imagesc(xx(1,:),yy(:,1),unormz(:,:)) ; shading interp
  axis square
  set(gca,'FontSize',14)
  xlabel('u_z','FontSize',14)
  colorbar
  subplot(2,2,4)
  imagesc(xx(1,:),yy(:,1),sqrt(unormx(:,:).^2+unormy(:,:).^2+unormz(:,:).^2)) ; shading interp
  axis square
  set(gca,'FontSize',14)
  xlabel('|u|','FontSize',14)
  colorbar
  figure('Name','2D (normal) log|u| slice')
  imagesc(xx(1,:),yy(:,1),log(sqrt(unormx(:,:).^2+unormy(:,:).^2+unormz(:,:).^2))) ; shading interp
  set(gca,'FontSize',14)
  xlabel('x','FontSize',14)
  ylabel('y','FontSize',14)
  colorbar
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n=s;
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
cutoff=10;
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
if std(sqrt(unormx.^2+unormy.^2+unormz.^2))>0
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
    dummy_spect=k.^(fit);
    scaling_factor=sum(spectnorm(1:10))/sum(dummy_spect(1:10));
    dummy_spect=dummy_spect*scaling_factor;
    hold on
    loglog(k,dummy_spect,'--k','LineWidth',2)
  end
end
