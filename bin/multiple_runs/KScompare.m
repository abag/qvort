function KScompare(runcount)
for i=1:runcount
  fIN1=sprintf('run%1d/data/KSwavenumbers.log',i);
  fIN2=sprintf('run%1d/data/ts.log',i);
  a=load(fIN1);
  rey(i)=(a(length(a))/a(1))^(4/3);
  b=load(fIN2);
  ll(i)=mean(b(int32(0.9*length(b)):length(b),6));
  %plot(b(:,6)) ; hold on
end
plot(rey,ll,'o')
