%a script to plot PDFs of loop sizes
function loop_size_pdf(start,finish,skip)
if nargin<3
  skip=1;
end
ts=load('./data/ts.log');
t=ts(:,2);
store_caxis=([min(t) max(t)]);
cmap=colormap(jet(finish)) ;
for i=start:skip:finish
  filename=sprintf('./data/loop_size%03d.log',i);
  A=load(filename);
  %[n xout]=ksdensity(A(:,2),'support','positive');
  [n xout]=ksdensity(A(:,2));
  plot(xout,n,'Color',cmap(i,:))
  hold on
end
hold off
caxis(store_caxis)
colorbar
set(gca,'FontSize',14)
ylabel('PDF','FontSize',14)
xlabel('N (points)','FontSize',14)
