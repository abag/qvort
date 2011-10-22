function KWC_amp(filenumber)
printit=1;
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
%%%%%%%%%%%%ORIGINAL DATA%%%%%%%%%%
subplot(5,1,1);
plot(amp(:,1),sqrt(amp(:,2).^2+amp(:,3).^2),'o');
title('original points')
xlabel('z','Fontsize',14)
ylabel('a','Fontsize',14)
set(gca,'Fontsize',14)
%%%%%%%%%%%%%SORT DATA%%%%%%%%%%%%%
subplot(5,1,2);
amp2=sortrows(amp,1);
plot(amp2(:,1),sqrt(amp2(:,2).^2+amp2(:,3).^2),'-','LineWidth',2);
title('ordered points')
xlabel('z','Fontsize',14)
ylabel('a','Fontsize',14)
set(gca,'Fontsize',14)
%%%%%%%%%%%%%INTERPOLATE DATA%%%%%%%
X=(-dims(2)/2:dims(1)/2:dims(2)/2.);
Yx=interp1(amp2(:,1),amp2(:,2),X);
Yy=interp1(amp2(:,1),amp2(:,3),X);
Z=Yx+i.*Yy;
subplot(5,1,3);
plot(X,abs(Z),'-','LineWidth',2);
title('interpolated to uniform mesh')
xlabel('z','Fontsize',14)
ylabel('a','Fontsize',14)
set(gca,'Fontsize',14)
%%%%%%%%%%%%%%%%PLOT PHASE%%%%%%%%%%%%%
subplot(5,1,4);
plot(X,angle(Z),'-','LineWidth',2);
title('phase')
xlabel('z','Fontsize',14)
ylabel('\theta','Fontsize',14)
set(gca,'Fontsize',14)
%%%%%%%%%%%%GRADIENT%%%%%%%%%%%%%
subplot(5,1,5);
plot(X,gradient(abs(Z),dims(2)/2),'-','LineWidth',2);
title('derivative')
xlabel('z','Fontsize',14)
ylabel('a dash','Fontsize',14)
set(gca,'Fontsize',14)
if printit==1
    print -depsc amp_info1.eps
end
%%%%%%%%%%%%%%%%%%REMOVE END POINTS%%%%%%%%%%%%%%%%%%
figure
X2=X(3:length(X)-3);
Z2=Z(3:length(Z)-3);
plot(X2,abs(Z2),'-','LineWidth',2);
xlabel('z','Fontsize',14)
ylabel('a','Fontsize',14)
set(gca,'Fontsize',14)
if printit==1
    print -depsc amp_info2.eps
end
%%%%%%%%%%%%%%SPECTRA%%%%%%%%%%%%%
N=length(X2);
P = abs(fft(Z2))/(N/2);
P = P(1:N/2).^2;
k=linspace(dims(1),4*pi/dims(2),N/2);
figure
loglog(k,P)
%------------compare with pwelch---------
figure
Pw=pwelch(Z2);
Pw=Pw(1:length(Pw)/2);
kw=linspace(dims(1),4*pi/dims(2),length(Pw));
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

