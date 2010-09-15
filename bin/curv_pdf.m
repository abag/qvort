%a script to plot histograms of curvature
%note I assume that the number of bins is 10
%run the script with the command 'print' to print to file
function curv_hist(option)
if nargin==0     
  option='empty';
end
switch option
case 'print'
    disp('will not print to screen but instead to .eps files')
case 'empty'
    otherwise
    disp('incorrect option, aborting script and printing help:')
    help curv_hist
    return
end
A=load('./data/curv_pdf.log');
ts=load('./data/ts.log');
t=ts(:,2);
s=size(A) ; snap_number=s(1)/10 ;
B=reshape(A,10,snap_number,2) ;
store_caxis=([min(t) max(t)]);
cmap=colormap(jet(snap_number)) ;
switch option
  case 'print'
    figure('visible','off');
  otherwise
    figure('Name', 'curvature PDF')      
end
for i=1:8:snap_number
 plot(B(:,i,1),B(:,i,2),'-','Color',cmap(i,:)) ;
 hold on   
end
hold off
xlabel('\kappa')
ylabel('PDF(\kappa)')
caxis(store_caxis)
colorbar
set(gca,'FontSize',14)
switch option
  case 'print'
    disp('printing to curv_pdf.eps')
    print('-depsc','./curv_pdf.eps')
end