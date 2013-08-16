function vel_bags(start,finish)
vel_stat_process(start,finish)
load ./velocity.mat   
ux2=abs(uy);
index=find(ux2<0.000001);
ux2(index)=[];
clear index
index=find(ux2>1.);
ux2(index)=[];
histnorm(ux2,100);
pause
[n xout]=histnorm(ux2,100);
plot(log(xout(50:end)),log(n(50:end)))
polyfit(log(xout(50:end)),log(n(50:end)),1)
