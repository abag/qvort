function acc_bags(start,finish,fit)
global xout
if nargin<3
    fit=0;
end
fit=1
delete('velocity.mat')
acc_stat_process(start,finish)
load velocity.mat
index=find(uy>0);
ux2=uy(index);
ux2=ux2/std(ux2);
[n xout]=histnorm(ux2,100);
pd=fitdist(ux2','Lognormal')
pdfA=lognpdf(xout,pd.mu,pd.sigma);
p=polyfit(log(xout(20:end)),log(n(20:end)),1);
if fit==1
  semilogy(xout,n,'ko',xout,pdfA,'b-',xout,2*exp(p(2))*xout.^-1.9,'r--','LineWidth',2)
else
  semilogy(xout,n,'ko','LineWidth',2)
end
xlabel('aa/aa','FontSize',16)
ylabel('PDF','FontSize',16)
set(gca,'FontSize',16)
%figure
%plot(log(xout(50:end)),log(n(50:end)))
polyfit(log(xout(50:end)),log(n(50:end)),1)
%pd=fitdist(ux2','Exponential')
