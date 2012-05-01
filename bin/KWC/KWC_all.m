function KWC_all(start,finish,skip)
counter=1;
k=0.;
P=0.;
kw=0.;
Pw=0.;
for i=start:skip:finish
  [dum_k dum_P dum_kw dum_Pw]=KWC_amp(i,0);
  k=k+dum_k;
  P=P+dum_P;
  kw=kw+dum_kw;
  Pw=Pw+dum_Pw;
  counter=counter+1;
end
counter=counter-1;
k=k/counter;
P=P'/counter;
kw=kw'/counter;
Pw=Pw/counter;
%%%%%%%%%%%%%%%%%%%%%%DO OUR SPECTRA%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
loglog(k,P,'LineWidth',2)
xlabel('log_{10} k','Fontsize',16)
ylabel('log_{10} A','Fontsize',16)
set(gca,'Fontsize',16)
length(k)
length(P)
p=polyfit(log10(k(30:100)),log10(P(30:100)),1)
figure('Name','Spectra with fit')
plot(log10(k(1:length(k)-100)),log10(P(1:length(k)-100)),'k','LineWidth',2)
hold on
plot(log10(k(15:length(k)-40)),p(1)*log10(k(15:length(k)-40))+p(2),'r--','LineWidth',1)
hold off
xlabel('log_{10} k','Fontsize',16)
ylabel('log_{10} A','Fontsize',16)
set(gca,'Fontsize',16)
figure('Name','Compensated -17/5')
plot(log(k(14:length(k)-100)),smooth(log((k(14:length(k)-100).^(17/5)).*P(14:length(k)-100)),10),'k','LineWidth',2)
figure('Name','Compensated -11/3')
plot(log(k(14:length(k)-100)),smooth(log((k(14:length(k)-100).^(11/3)).*P(14:length(k)-100)),10),'k','LineWidth',2)
figure('Name','Componstated on same scale')
plot(log10(k(10:length(k)-100)),smooth(log10((k(10:length(k)-100).^(17/5)).*P(10:length(k)-100)),10),'k','LineWidth',2)
hold on
plot(log10(k(10:length(k)-100)),smooth(log10((k(10:length(k)-100).^(11/3)).*P(10:length(k)-100)),10),'r','LineWidth',2)
hold on
plot(log10(k(10:length(k)-100)),smooth(log10((k(10:length(k)-100).^(3)).*P(10:length(k)-100)),10),'b','LineWidth',2)
xlabel('log_{10} k','Fontsize',16)
ylabel('log_{10} Ak^3, Ak^{3.4}, Ak^{3.66} ','Fontsize',16)
set(gca,'Fontsize',16)
return
%%%%%%%%%%%%%%%%%%%%%%NOW WELCH%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p2=polyfit(log(kw(14:length(kw)-60)),log(Pw(14:length(kw)-60)),1)
figure('Name','Welch spectra with fit')
plot(log(kw(14:length(kw)-20)),log(Pw(14:length(kw)-20)),'k','LineWidth',2)
hold on
plot(log(kw(14:length(kw)-20)),p2(1)*log(kw(14:length(kw)-20))+p2(2),'r--','LineWidth',1)
plot(log(kw(14:length(kw)-20)),-3.66.*log(kw(14:length(kw)-20))+p2(2)-1.5,'b--','LineWidth',1)
plot(log(kw(14:length(kw)-20)),-3.4.*log(kw(14:length(kw)-20))+p2(2)-2,'m--','LineWidth',1)
xlabel('log k','Fontsize',16)
ylabel('log A','Fontsize',16)
figure('Name','Welch Compensated -17/5')
plot(log(kw(10:length(kw)-50)),log((kw(10:length(kw)-50).^(17/5)).*Pw(10:length(kw)-50)),'k','LineWidth',2)
figure('Name','Welch Compensated -11/3')
plot(log(kw(10:length(kw)-50)),log((kw(10:length(kw)-50).^(11/3)).*Pw(10:length(kw)-50)),'k','LineWidth',2)
hold on
plot(log(kw(9:length(kw)-40)),log(0.0003*(kw(9:length(kw)-40).^(11/3)).*(kw(9:length(kw)-40).^(-3.4))),'r--','LineWidth',1)

