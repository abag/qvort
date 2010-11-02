%a function which loops over all available var files and computes the correlation
%dimension of the vortex filaments, then plots this as a function of time.
function coor_dim_all(skip)
  if nargin==0
    skip=1;
  end 
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
      p=corr_dim(i,'noplot');
      slope(i)=p(1);
  end
  plot(t(1:skip:s),slope(1:skip:s),'k','LineWidth',2)
  set(gca,'FontSize',14)
  xlabel('t','FontSize',14)
  ylabel('corr','FontSize',14)
