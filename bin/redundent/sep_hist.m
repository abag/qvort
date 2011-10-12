%a script to plot histograms of particle setparation
%run the script with the command 'print' to print to file
function sep_hist(start,final,skip,option)
if nargin==3     
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch option
  case 'print'
    figure('visible','off');
  otherwise
    figure('Name', 'separation PDF')      
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=start:skip:final
  [f xi]=vortex_sep_hist(i,'noplot');
  switch option
    case 'loglog'
      loglog(B(:,i,1),log(B(:,i,2)),'-','Color',cmap(i,:)) ;
    case 'log'
      semilogx(xi,f,'-','Color',cmap(i,:)) ; 
    otherwise
      plot(xi,f,'-','Color',cmap(i,:)) ; 
  end
  hold on   
end
hold off
switch option
  case 'loglog'
    xlabel('log sep','FontSize',14)
    ylabel('log PDF(sep)','FontSize',14)
  case 'log'
    xlabel('sep','FontSize',14)
    ylabel('log PDF(sep)','FontSize',14)
  otherwise
    xlabel('sep','FontSize',14)
    ylabel('PDF(sep)','FontSize',14)
end
caxis(store_caxis)
colorbar
set(gca,'FontSize',14)
switch option
  case 'print'
    disp('printing to sep_pdf.eps')
    print('-depsc','./sep_pdf.eps')
end
