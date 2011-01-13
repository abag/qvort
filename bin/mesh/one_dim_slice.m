function one_dim_slice(filenumber,fit,fitnorm)
if nargin<2
  do_fit=0;
else
  do_fit=1;
end
if nargin<3
  fitnorm=fit;
end
filename=sprintf('data/vel_slice_1D%03d.log',filenumber);
A=load(filename);
u2=sqrt(A(:,2).^2+A(:,3).^2+A(:,4).^2);
if std(u2)>0.
  figure('Name','1D (superfluid) velocity slices')
  subplot(4,1,1)
    plot(A(:,1),A(:,2),'LineWidth',2,'Color','m')
    set(gca,'FontSize',14)
    ylabel('u_x','FontSize',14)
  subplot(4,1,2)
    plot(A(:,1),A(:,3),'LineWidth',2,'Color','b')
    set(gca,'FontSize',14)
    ylabel('u_y','FontSize',14)
  subplot(4,1,3)
    plot(A(:,1),A(:,4),'LineWidth',2,'Color','r')
    set(gca,'FontSize',14)
    ylabel('u_z','FontSize',14)
    xlabel('x','FontSize',14)
  subplot(4,1,4)
    plot(A(:,1),u2,'LineWidth',2,'Color','k')
    set(gca,'FontSize',14)
    ylabel('|u|','FontSize',14)
    xlabel('x','FontSize',14)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u2norm=sqrt(A(:,5).^2+A(:,6).^2+A(:,7).^2);
if std(u2norm)>0.
  figure('Name','1D (normal) velocity slices')
  subplot(3,1,1)
    plot(A(:,1),A(:,5),'LineWidth',2,'Color','g')
    set(gca,'FontSize',14)
    ylabel('u_{norm,x}','FontSize',14)
  subplot(3,1,2)
    plot(A(:,1),A(:,6),'LineWidth',2,'Color','k')
    set(gca,'FontSize',14)
    ylabel('u_{norm,y}','FontSize',14)
  subplot(3,1,3)
    plot(A(:,1),A(:,7),'LineWidth',2,'Color','b')
    set(gca,'FontSize',14)
    ylabel('u_{norm,z}','FontSize',14)
    xlabel('x','FontSize',14)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Name','1D (superfluid) velocity spectrum')
n=length(A);
fux=fftn(A(:,2))/n;
fuy=fftn(A(:,3))/n;
fuz=fftn(A(:,4))/n;
energyr=real(fux).^2+real(fuy).^2+real(fuz).^2;
energyi=imag(fux).^2+imag(fuy).^2+imag(fuz).^2;
midpt=n/2+1;
spect(1:1.5*n)=0.;
for i=1:n
  ii=i;
  if ii>midpt ; ii=n-ii+1; ; end ;
  spect(ii)=spect(ii)+energyr(i)+energyi(i);
end
k=1:midpt-1;
loglog(k,spect(1:midpt-1),'LineWidth',2)
  set(gca,'FontSize',14)
  ylabel('log E(k)','FontSize',14)
  xlabel('log k','FontSize',14)
if do_fit==1
  disp(sprintf('fitting a slope of %f to plot',fit))
  dummy_spect=k.^(fit);
  scaling_factor=sum(spect(1:10))/sum(dummy_spect(1:10));
  dummy_spect=dummy_spect*scaling_factor;
  hold on
  loglog(k,dummy_spect,'--k','LineWidth',2)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if std(u2norm)>0.
  figure('Name','1D (normal) velocity spectrum')
  fux=fftn(A(:,5))/n;
  fuy=fftn(A(:,6))/n;
  fuz=fftn(A(:,7))/n;
  energyr=real(fux).^2+real(fuy).^2+real(fuz).^2;
  energyi=imag(fux).^2+imag(fuy).^2+imag(fuz).^2;
  spectnorm(1:1.5*n)=0.;
  for i=1:n
    ii=i;
    if ii>midpt ; ii=n-ii+1; ; end ;
    spectnorm(ii)=spectnorm(ii)+energyr(i)+energyi(i);
  end
  k=1:midpt-1;
  loglog(k,spectnorm(1:midpt-1),'-r','LineWidth',2)
  set(gca,'FontSize',14)
  ylabel('log E(k)','FontSize',14)
  xlabel('log k','FontSize',14)
  if do_fit==1
    disp(sprintf('fitting a slope of %f to plot',fitnorm))
    dummy_spect=k.^(fitnorm);
    scaling_factor=sum(spectnorm(1:10))/sum(dummy_spect(1:10));
    dummy_spect=dummy_spect*scaling_factor;
    hold on
    loglog(k,dummy_spect,'--k','LineWidth',2)
  end
end
