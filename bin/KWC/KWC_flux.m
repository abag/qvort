function KWC_all(start,finish)
A=load('./data/ts.log');
counter=1;
k=0.;
P=0.;
kw=0.;
Pw=0.;
for i=start:finish
  [dum_k dum_P dum_kw dum_Pw]=KWC_amp(i,0);
  k=dum_k;
  P(counter,1:length(k))=dum_P;
  counter=counter+1;
end
dims=load('./data/dims.log');
k=k*2*pi/dims(1);
cmap=colormap(jet(size(P,2)));
figure('Name','Time series n_k(t)')
for i=1:30:size(P,2)
  semilogy(A(start:finish,2),P(:,i),'Color',cmap(i,:))
  hold on
end
hold off
figure('Name','Time series dn_k/dt')
for i=1:size(P,2)
Pdot(:,i)=gradient(P(:,i),A(2,2)-A(1,2));
end
for i=1:30:size(P,2)
  plot(A(start:finish,2),Pdot(:,i),'Color',cmap(i,:))
  hold on
end
hold off
figure('Name','Time series int dn_k/dt dk')
for i=1:size(P,1)
intP(i,1)=0.;
i
for j=2:size(P,2)
  intP(i,j)=trapz(Pdot(i,1:j)',k(1:j));
end
end
mean_intP=mean(intP,1);
semilogx(mean_intP)


