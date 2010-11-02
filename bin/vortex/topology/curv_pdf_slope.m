%a script to plot histograms of curvature
%note I assume that the number of bins is 10
%run the script with the command 'print' to print to file
function curv_hist
A=load('./data/curv_pdf.log');
ts=load('./data/ts.log');
t=ts(:,2);
s=size(A) ; snap_number=s(1)/10 ;
B=reshape(A,10,snap_number,2) ;
figure('Name', 'curvature PDF slope')      
for i=1:snap_number
  if B(6,i,2)==0
    p=polyfit(log(B(2:5,i,1)),log(B(2:5,i,2)),1);
  elseif B(7,i,2)==0
    p=polyfit(log(B(3:6,i,1)),log(B(3:6,i,2)),1);
  elseif B(8,i,2)==0
    p=polyfit(log(B(4:7,i,1)),log(B(4:7,i,2)),1);
  else
    p=polyfit(log(B(4:8,i,1)),log(B(4:8,i,2)),1);
  end
  curv_slope(i)=p(1);
end
set(gca,'FontSize',14)
 s=length(t);
 c1=smooth(curv_slope,100)';
 c2=smooth(curv_slope,1000)';
 plot(t(1:s-200),c1(1:s-200),'k')
 hold on
 plot(t(1:s-200),c2(1:s-200),'--r','LineWidth',2)
 xlabel('t','FontSize',14)
 ylabel('slope','FontSize',14)
 mc=mean(curv_slope(floor(0.1*s):s-1));
 disp(sprintf('mean of the curvature pdf slope: %06f',mc));