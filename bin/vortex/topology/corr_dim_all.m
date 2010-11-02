function coor_dim_all
  ts=load('./data/ts.log');
  t=ts(:,2);
  s=length(t);
  for i=1:s
      p=corr_dim(i,'noplot');
      slope(i)=p(1);
  end
  plot(t,slope,'k','LineWidth',2)
  set(gca,'FontSize',14)
  xlabel('t','FontSize',14)
  ylabel('corr','FontSize',14)
