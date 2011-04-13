function one_dim_slice_avg(start,final,skip)
if nargin<1
  start=1 ;
end
if nargin<2
  final=0;
end
if nargin<3
  skip=1;
end
A=load('data/ts.log');
filenumber=length(A);
if (final==0) || (final>filenumber)
  final=filenumber;
end
fprintf('starting file: %d final file: %d \n',start,final)
if (skip>1)
  fprintf('only processing every %d file \n',skip)
end
counter=1;
for i=start:skip:final
  filename=sprintf('data/vel_slice_1D%03d.log',i);
  A=load(filename);
  u2=sqrt(A(:,2).^2+A(:,3).^2+A(:,4).^2);
  n=length(A);
  fux=fftn(A(:,2))/n;
  fuy=fftn(A(:,3))/n;
  fuz=fftn(A(:,4))/n;
  energyr=real(fux).^2+real(fuy).^2+real(fuz).^2;
  energyi=imag(fux).^2+imag(fuy).^2+imag(fuz).^2;
  midpt=n/2+1;
  if i==start
    spect(1:1.5*n)=0.;
  end
  for i=1:n
    ii=i;
    if ii>midpt ; ii=n-ii+1; ; end ;
    spect(ii)=spect(ii)+energyr(i)+energyi(i);
  end
  counter=counter+1;
  clear A
end
spect=spect/(counter-1);
k=1:midpt-4;
figure('Name','Super')
loglog(k,spect(1:midpt-4),'LineWidth',2)
set(gca,'FontSize',14)
ylabel('log E(k)','FontSize',14)
xlabel('log k','FontSize',14)
%FIT TO UPPER END OF SPECTRA
  fit=-5/3;
  disp(sprintf('fitting a slope of %f to plot',fit))
  dummy_spect=k.^(fit);
  scaling_factor=sum(spect(1:10))/sum(dummy_spect(1:10));
  dummy_spect=dummy_spect*scaling_factor;
  hold on
  loglog(k,dummy_spect,'--k','LineWidth',2)
%FIT TO LOWER END
  fit=-1;
  disp(sprintf('fitting a slope of %f to plot',fit))
  dummy_spect=k.^(fit);
  scaling_factor=sum(spect(30:50))/sum(dummy_spect(30:50));
  dummy_spect=dummy_spect*scaling_factor;
  hold on
  loglog(k,dummy_spect,'--r','LineWidth',2)
%NOW CLEAR ALL AND DO NORMAL FLUID
counter=1;
for i=start:skip:final
  filename=sprintf('data/vel_slice_1D%03d.log',i);
  A=load(filename);
  u2=sqrt(A(:,5).^2+A(:,6).^2+A(:,7).^2);
  if u2<=0
    disp('no normal fluid, exiting')
    return
  end
  n=length(A);
  fux=fftn(A(:,5))/n;
  fuy=fftn(A(:,6))/n;
  fuz=fftn(A(:,7))/n;
  energyr=real(fux).^2+real(fuy).^2+real(fuz).^2;
  energyi=imag(fux).^2+imag(fuy).^2+imag(fuz).^2;
  midpt=n/2+1;
  if i==start
    spect(1:1.5*n)=0.;
  end
  for i=1:n
    ii=i;
    if ii>midpt ; ii=n-ii+1; ; end ;
    spect(ii)=spect(ii)+energyr(i)+energyi(i);
  end
  clear A
  counter=counter+1;
end
spect=spect/(counter-1);
k=1:midpt-4;
figure('Name','Normal')
loglog(k,spect(1:midpt-4),'LineWidth',2)
set(gca,'FontSize',14)
ylabel('log E(k)','FontSize',14)
xlabel('log k','FontSize',14)
%FIT TO UPPER END OF SPECTRA
  fit=-5/3;
  disp(sprintf('fitting a slope of %f to plot',fit))
  dummy_spect=k.^(fit);
  scaling_factor=sum(spect(1:10))/sum(dummy_spect(1:10));
  dummy_spect=dummy_spect*scaling_factor;
  hold on
  loglog(k,dummy_spect,'--k','LineWidth',2)
