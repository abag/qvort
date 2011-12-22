function [k P kw Pw]=KWC_amp(filenumber,plotit);
printit=0 ; 
%get the dimensions information from dims.log
dims=load('./data/dims.log');
filename=sprintf('data/var%04d.log',filenumber);
if dims(4)==1
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
  fclose(fid);
else 
  fid=fopen(filename);
  if fid<0
      disp('var file does not exist, exiting script')
      return
  end
  %read the time
  tline=fgetl(fid);
  dummy=textscan(tline, '%f');
  time=dummy{:};
  %how many particles
  tline=fgetl(fid);
  dummy=textscan(tline, '%d');  
  number_of_particles=dummy{:};
  %get the particles%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  for j=1:number_of_particles
    tline=fgetl(fid);
    dummy=textscan(tline, '%f');
    dummy_vect=dummy{:};
    x(j)=dummy_vect(1);
    y(j)=dummy_vect(2);
    z(j)=dummy_vect(3);
    f(j)=dummy_vect(4);
    u(j)=dummy_vect(5);
    u2(j)=dummy_vect(6);
  end
  f=uint16(f);
  fclose(fid);
end
counter=1;
for j=1:number_of_particles
  if f(j)>0
    amp(counter,1)=z(j);
    amp(counter,2)=x(j);
    amp(counter,3)=y(j);
    counter=counter+1;
  end
end
amp2=sortrows(amp,1);
X=(-dims(2)/2:dims(1)/2:dims(2)/2.);
Yx=interp1(amp2(:,1),amp2(:,2),X,'cubic');
Yy=interp1(amp2(:,1),amp2(:,3),X,'cubic');
Z=Yx+i.*Yy;
X2=X(5:length(X)-3);
Z2=Z(5:length(Z)-3);
if plotit==1
figure
plot(amp(:,1),sqrt(amp(:,2).^2+amp(:,3).^2),'o');
hold on
plot(X,abs(Z),'-','LineWidth',2);
figure
%%%%%%%%%%%%ORIGINAL DATA%%%%%%%%%%
fonts=10
subplot(5,1,1);
plot(amp(:,1),sqrt(amp(:,2).^2+amp(:,3).^2),'o');
title('original points')
xlabel('z','Fontsize',fonts)
ylabel('a','Fontsize',fonts)
set(gca,'Fontsize',fonts)
%%%%%%%%%%%%%SORT DATA%%%%%%%%%%%%%
subplot(5,1,2);
plot(amp2(:,1),sqrt(amp2(:,2).^2+amp2(:,3).^2),'-','LineWidth',2);
title('ordered points')
xlabel('z','Fontsize',fonts)
ylabel('a','Fontsize',fonts)
set(gca,'Fontsize',fonts)
%%%%%%%%%%%%%INTERPOLATE DATA%%%%%%%
subplot(5,1,3);
plot(X,abs(Z),'-','LineWidth',2);
title('interpolated to uniform mesh')
xlabel('z','Fontsize',fonts)
ylabel('a','Fontsize',fonts)
set(gca,'Fontsize',fonts)
%%%%%%%%%%%%%%%%PLOT PHASE%%%%%%%%%%%%%
subplot(5,1,4);
plot(X,angle(Z),'-','LineWidth',2);
title('phase')
xlabel('z','Fontsize',fonts)
ylabel('\theta','Fontsize',fonts)
set(gca,'Fontsize',fonts)
%%%%%%%%%%%%GRADIENT%%%%%%%%%%%%%
subplot(5,1,5);
plot(X,gradient(abs(Z),dims(1)/2),'-','LineWidth',2);
disp('gradient mean')
mean(abs(gradient(abs(Z),dims(1)/2)))
title('derivative')
xlabel('z','Fontsize',fonts)
ylabel('a dash','Fontsize',fonts)
set(gca,'Fontsize',fonts)
if printit==1
    print -depsc amp_info1.eps
end
%%%%%%%%%%%%%%%%%%REMOVE END POINTS%%%%%%%%%%%%%%%%%%
figure
plot(X2,abs(Z2),'-','LineWidth',2);
xlabel('z','Fontsize',14)
ylabel('a','Fontsize',14)
set(gca,'Fontsize',14)
if printit==1
    print -depsc amp_info2.eps
end
end
%%%%%%%%%%%%%%SPECTRA%%%%%%%%%%%%%
N=length(X2);
P = abs(fft(Z2))/(N/2);
P = P(1:N/2).^2;
%k=linspace(dims(1),4*pi/dims(2),N/2);
k=linspace(1,N/2,N/2);
Pw=pwelch(Z2);
Pw=Pw(1:length(Pw)/2);
kw=linspace(dims(1),4*pi/dims(2),length(Pw));
if plotit==1
%------------compare with pwelch---------
figure
loglog(k,P)
figure
loglog(Pw,'b')
hold on
loglog(P,'g')
if printit==1
    print -depsc spectra1.eps
end
%------------now work out scaling---------
figure
plot(log(k(20:length(k))),log(P(20:length(k))),'LineWidth',1.5)
p=polyfit(log(k(20:length(k))),log(P(20:length(k))),1)
hold on
fit1=p(1).*log(k(20:length(k)))+p(2);
fit2=-3.66.*log(k(20:length(k)))+p(2);
plot(log(k(20:length(k))),fit1,'k')
plot(log(k(20:length(k))),fit2,'r')
xlabel('log k','Fontsize',14)
ylabel('log A','Fontsize',14)
set(gca,'Fontsize',14)
if printit==1
    print -depsc spectra_fit.eps
end
figure
plot(log(kw(20:length(kw))),log(Pw(20:length(kw))),'LineWidth',2)
xlabel('log k','Fontsize',14)
ylabel('log A','Fontsize',14)
set(gca,'Fontsize',14)
end
