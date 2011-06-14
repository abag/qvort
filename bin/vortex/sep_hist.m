%a script to plot histograms of particle setparation
%run the script with the command 'print' to print to file
function sep_hist(option)
if nargin==0     
  option='empty';
end
switch option
case 'loglog'
    disp('will plot both axis on a logscale')
case 'log'
    disp('will plot with a logscale on the y-axis')
case 'print'
    disp('will not print to screen but instead to .eps files')
case 'empty'
    otherwise
    disp('incorrect option, aborting script and printing help:')
    help curv_hist
    return
end
ts=load('./data/ts.log');
t=ts(:,2);
s=size(ts);
store_caxis=([min(t) max(t)]);
cmap=colormap(jet(s(1))) ;
for i=1:s(1)
  [f xi]=vortex_sep_hist(i);
  plot(xi,f,'Color',cmap(i,:)) ; 
end 
%switch option
%  case 'print'
%    figure('visible','off');
%  otherwise
%    figure('Name', 'curvature PDF')      
%end
%for i=1:8:snap_number
%  switch option
%    case 'loglog'
%      loglog(B(:,i,1),B(:,i,2),'-','Color',cmap(i,:)) ;
%      xlabel('log \kappa','FontSize',14)
%      ylabel('log PDF(\kappa)','FontSize',14)
%    case 'log'
%      plot(B(:,i,1),log(B(:,i,2)),'-','Color',cmap(i,:)) ;
%      xlabel('\kappa','FontSize',14)
%      ylabel('log PDF(\kappa)','FontSize',14)
%    otherwise
%      plot(B(:,i,1),B(:,i,2),'-','Color',cmap(i,:)) ;
%      %plot(B(:,i,1),B(:,i,2),'-','LineWidth',2,'Color',cmap(i,:)) ;
%      xlabel('\kappa','FontSize',14)
%      ylabel('PDF(\kappa)','FontSize',14)
%  end
%  hold on   
%end
%hold off
%caxis(store_caxis)
%colorbar
%set(gca,'FontSize',14)
%switch option
%  case 'print'
%    disp('printing to curv_pdf.eps')
%    print('-depsc','./curv_pdf.eps')
%end
