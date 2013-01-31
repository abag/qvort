function KWC_all(start,finish,skip)
cmap=colormap(jet(finish));
for i=start:skip:finish
  [dum_k dum_P dum_kw dum_Pw]=KWC_amp(i,0);
  loglog(dum_k,dum_P,'LineWidth',2,'Color',cmap(i,:))
  hold on
end
set(gca,'Fontsize',16)


