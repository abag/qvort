%a function which loops over all available var files and computes the fractal
%dimension of the vortex filaments, then plots this as a function of time.
%takes the option method which should be set to haus or corr
function coor_dim_all(start,finish,skip,method)
  if nargin<3
    disp('need thre input arguments start,finish,skip')
    return
  end 
  if nargin==3
    method='haus';
  end 
  ts=load('./data/ts.log');
  t=ts(:,2);
  s=length(t);
  for i=start:skip:finish
      if skip<10
        if mod(i,10)==0
          disp(sprintf('processing varfile %04d of %04d',i,s))
        end
      else
        disp(sprintf('processing varfile %04d of %04d',i,s))
      end
      switch method
        case 'haus'
      p=haus_dim(i);
      slope(i)=p;
        case 'corr'
        p=corr_dim(i,'noplot');
        slope(i)=p(1);
      end
  end
  plot(t(start:skip:finish),slope(start:skip:finish),'k','LineWidth',2)
  set(gca,'FontSize',14)
  xlabel('t','FontSize',14)
  ylabel('corr','FontSize',14)
  t2=t(start:skip:finish) ; correlation_dim=slope(start:skip:finish);
  save correlation_dimension.mat t2 correlation_dim
