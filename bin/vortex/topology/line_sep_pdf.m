%a script to plot histograms of point separations
%note I assume that the number of bins is 20
%run the script with the command 'print' to print to file
function line_sep_hist(option)
%set time limits to plot
min_t=1.95 ; max_t=2.2
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
A=load('./data/line_sep_pdf.log');
ts=load('./data/ts.log');
t=ts(:,2);
bin_number=10;
s=size(A) ; snap_number=s(1)/bin_number ;
B=reshape(A,bin_number,snap_number,2) ;
store_caxis=([min(t) max(t)]);
t_plot=max(t)/snap_number;
cmap=colormap(jet(snap_number)) ;
switch option
  case 'print'
    figure('visible','off');
  otherwise
    figure('Name', 'separation PDF')      
end
for i=1:snap_number
  if (i*t_plot)>min_t && (i*t_plot)<max_t
   %pause and give time information below
   %i*t_plot
   %pause
    switch option
      case 'loglog'
        loglog(B(:,i,1),B(:,i,2),'-','Color',cmap(i,:)) ;
        xlabel('log line sep','FontSize',14)
        ylabel('log PDF(line sep)','FontSize',14)
      case 'log'
        semilogy(B(:,i,1),B(:,i,2),'-','Color',cmap(i,:)) ;
        xlabel('line sep','FontSize',14)
        ylabel('log PDF(line sep)','FontSize',14)
      otherwise
        plot(B(:,i,1),B(:,i,2),'-','Color',cmap(i,:)) ;
        %thicker lines
        %plot(B(:,i,1),B(:,i,2),'-','LineWidth',2,'Color',cmap(i,:)) ;
        xlabel('line sep','FontSize',14)
        ylabel('PDF(line sep)','FontSize',14)
    end
  end
  hold on   
end
hold off
caxis(store_caxis)
colorbar
set(gca,'FontSize',14)
switch option
  case 'print'
    disp('printing to line_sep_pdf.eps')
    print('-depsc','./line_sep_pdf.eps')
end
