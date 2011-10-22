function KWC_all(start,finish,skip)
counter=1
k=0.;
P=0.;
for i=start:skip:finish
  [dummy_k dummy_P]=KWC_amp(i);
  k=k+dummy_k;
  P=P+dummy_P;
  counter=counter+1
end
k=k/counter;
P=P/counter;
loglog(k,P)
