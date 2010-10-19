function KScompare(runcount)
for i=1:runcount %get reynolds numbers first
  fIN1=sprintf('run%1d/data/KSwavenumbers.log',i);
  a=load(fIN1);
  rey(i)=(a(length(a))/a(1))^(4/3);
end
store_caxis=([min(rey) max(rey)]);
range=ceil(max(rey)-min(rey))+1;
cmap=colormap(jet(range)) ;
for i=1:runcount
    fIN2=sprintf('run%1d/data/ts.log',i);
    b=load(fIN2);
    subplot(1,2,1)
      plot(b(:,6),'LineWidth',2,'Color',cmap(ceil(rey(i)-min(rey)+1),:)) ; hold on
     ll(i)=mean(b(int32(0.9*length(b)):length(b),6));
end
hold off
xlabel('t','FontSize',12)
ylabel('line length','FontSize',12)
%caxis(store_caxis)
%colorbar
axis tight
set(gca,'FontSize',14)
for i=1:runcount
    subplot(1,2,2)
    plot(rey(i),ll(i),'o','MarkerSize',6,'MarkerFaceColor',cmap(ceil(rey(i)-min(rey)+1),:),'MarkerEdgeColor','k')
    hold on
end
hold off
xlabel('Re','FontSize',12)
ylabel('line length','FontSize',12)
set(gca,'FontSize',14)
