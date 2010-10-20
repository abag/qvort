function vortex_smooth_slope
  ts=load('./data/ts.log');
  t=ts(:,2);
  s=length(t);
  for i=1:s
      p=vortex_smooth(i,10,'noplot');
      slope(i)=p(1);
  end
  plot(t,slope,'k','LineWidth',2)
  set(gca,'FontSize',14)
  xlabel('t','FontSize',14)
  ylabel('slope','FontSize',14)