function acc_bags(start,finish)
acc_stat_process(start,finish)
load ./velocity.mat   
index=find(ux>0);
ux2=ux(index);
ux2=ux2/std(ux2);
[n xout]=histnorm(ux2,200);
pd=fitdist(ux2','Lognormal')
pdfA=lognpdf(xout,pd.mu,pd.sigma);
p=polyfit(log(xout(50:end)),log(n(50:end)),1);

semilogy(xout,n,'o',xout,pdfA,xout,exp(p(2))*xout.^p(1))
%figure
%plot(log(xout(50:end)),log(n(50:end)))
polyfit(log(xout(50:end)),log(n(50:end)),1)
%pd=fitdist(ux2','Exponential')
