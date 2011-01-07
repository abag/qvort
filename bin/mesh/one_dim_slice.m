function one_dim_slice(fit)
if nargin<1
  do_fit=0
else
  do_fit=1
end
A=load('data/vel_slice_1D.log');
figure('Name','1D (superfluid) velocity slices')
subplot(3,1,1)
  plot(A(:,1),A(:,2),'LineWidth',2,'Color','m')
  set(gca,'FontSize',14)
  ylabel('u_x','FontSize',14)
subplot(3,1,2)
  plot(A(:,1),A(:,3),'LineWidth',2,'Color','b')
  set(gca,'FontSize',14)
  ylabel('u_y','FontSize',14)
subplot(3,1,3)
  plot(A(:,1),A(:,4),'LineWidth',2,'Color','r')
  set(gca,'FontSize',14)
  ylabel('u_z','FontSize',14)
  xlabel('x','FontSize',14)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
  scaling_factor=spect(2)/dummy_spect(2);
  dummy_spect=dummy_spect*scaling_factor;
  hold on
  loglog(k,dummy_spect,'--k','LineWidth',2)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
fitt=1
if fitt==1
  dummy_spect=k.^(-5/3);
  scaling_factor=spectnorm(2)/dummy_spect(2);
  dummy_spect=dummy_spect*scaling_factor;
  hold on
  loglog(k,dummy_spect,'--k','LineWidth',2)
end
