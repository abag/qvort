%a script to plot histograms of point separations
%note I assume that the number of bins is 20
%run the script with the command 'print' to print to file
function curv_hist(option)
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
A=load('./data/sep_pdf.log');
ts=load('./data/ts.log');
t=ts(:,2);
bin_number=20;
s=size(A) ; snap_number=s(1)/bin_number ;
B=reshape(A,bin_number,snap_number,2) ;
store_caxis=([min(t) max(t)]);
cmap=colormap(jet(snap_number)) ;
switch option
  case 'print'
    figure('visible','off');
  otherwise
    figure('Name', 'separation PDF')      
end
for i=1:snap_number
  switch option
    case 'loglog'
      loglog(B(:,i,1),B(:,i,2),'-','Color',cmap(i,:)) ;
      xlabel('log sep','FontSize',14)
      ylabel('log PDF(sep)','FontSize',14)
    case 'log'
      semilogy(B(:,i,1),B(:,i,2),'-','Color',cmap(i,:)) ;
      xlabel('sep','FontSize',14)
      ylabel('log PDF(sep)','FontSize',14)
    otherwise
      plot(B(:,i,1),B(:,i,2),'-','Color',cmap(i,:)) ;
      %plot(B(:,i,1),B(:,i,2),'-','LineWidth',2,'Color',cmap(i,:)) ;
      xlabel('sep','FontSize',14)
      ylabel('PDF(sep)','FontSize',14)
  end
  hold on   
end
hold off
caxis(store_caxis)
colorbar
set(gca,'FontSize',14)
switch option
  case 'print'
    disp('printing to sep_pdf.eps')
    print('-depsc','./sep_pdf.eps')
end
