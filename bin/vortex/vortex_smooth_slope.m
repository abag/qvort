%a function which loops over all available var files and computes the slope of the 
%ampitude spectrum, then plots this as a function of time.
%need to add in vortex angle call here
function vortex_smooth_slope(skip)
  ts=load('./data/ts.log');
  t=ts(:,2);
  s=length(t);
  for i=1:skip:s
      if skip<10
        if mod(i,10)==0
          disp(sprintf('processing varfile %04d of %04d',i,s))
        end
      else
        disp(sprintf('processing varfile %04d of %04d',i,s))
      end
      [p meana maxa]=vortex_smooth(i,10,'noplot');
      slope(i)=p(1);
      abar(i)=meana;
      amax(i)=maxa;
  end
  figure('Name','spectrum slope')
    plot(t(1:skip:s),slope(1:skip:s),'k','LineWidth',2)
    %plot(t(1:skip:s),smooth(slope(1:skip:s),3),'k','LineWidth',2)
    set(gca,'FontSize',14)
    xlabel('t','FontSize',14)
    ylabel('slope','FontSize',14)
  figure('Name','mean amplitude')
    plot(t(1:skip:s),abar(1:skip:s),'k','LineWidth',2)
    %plot(t(1:skip:s),smooth(abar(1:skip:s),3),'k','LineWidth',2)
    set(gca,'FontSize',14)
    xlabel('t','FontSize',14)
    ylabel('meana','FontSize',14)
end
