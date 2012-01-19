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
  P(counter,:)=dum_P;
  counter=counter+1;
end
semilogy(A(start:end,6),P)

