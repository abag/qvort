%a script to plot histograms of curvature
%note I assume that the number of bins is 10
%run the script with the command 'print' to print to file
%run the script with the command 'loglog' for double log scale
%run the script with the command 'log' for semi-log (y axis)
function curv_pdf(option)
if nargin==0     
  option='empty';
end
switch option
case 'loglog'
    disp('will plot both axis on a logscale')
case 'logy'
    disp('will plot with a logscale on the y-axis')
case 'logx'
    disp('will plot with a logscale on the x-axis')
case 'print'
    disp('will not print to screen but instead to .eps files')
case 'empty'
    otherwise
    disp('incorrect option, aborting script and printing help:')
    help curv_hist
    return
end
ignore_percent=0.5;
A=load('./data/curv_pdf.log');
ts=load('./data/ts.log');
t=ts(:,2);
s=size(A) ; snap_number=s(1)/10 ;
B=reshape(A,10,snap_number,2) ;

switch option
  case 'print'
    figure('visible','off');
  otherwise
    figure('Name', 'curvature PDF')      
end
avg_x=0.;
avg_y=0.;
floor(snap_number*ignore_percent)
for i=floor(snap_number*ignore_percent):snap_number
  avg_x=avg_x+B(:,i,1);
  avg_y=avg_y+B(:,i,2);
end
avg_x=avg_x/(snap_number-floor(snap_number*ignore_percent));
avg_y=avg_y/(snap_number-floor(snap_number*ignore_percent));
switch option
    case 'loglog'
      loglog(avg_x,avg_y,'LineWidth',2) ;
      xlabel('log \kappa','FontSize',14)
      ylabel('log PDF(\kappa)','FontSize',14)
    case 'logy'
      semilogy(avg_x,avg_y,'LineWidth',2) ;
      xlabel('\kappa','FontSize',14)
      ylabel('log PDF(\kappa)','FontSize',14)
    case 'logx'
      semilogx(avg_x,avg_y,'LineWidth',2) ;
      xlabel('log \kappa','FontSize',14)
      ylabel('PDF(\kappa)','FontSize',14)
    otherwise
      plot(avg_x,avg_y,'LineWidth',2) ;
      %plot(B(:,i,1),B(:,i,2),'-','LineWidth',2,'Color',cmap(i,:)) ;
      xlabel('\kappa','FontSize',14)
      ylabel('PDF(\kappa)','FontSize',14)
end
set(gca,'FontSize',14)
switch option
  case 'print'
    disp('printing to curv_pdf.eps')
    print('-depsc','./curv_pdf.eps')
end
save avg_curv_pdf.mat avg_x avg_y