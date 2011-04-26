function vortex_anim(start,final,skip)
if nargin<1
  disp('I at least need finish filenumbers')
  return
elseif nargin<2
  disp('Assuming number given is final, start set to 1')
  final=start;
  start=1;
  skip=1.
elseif nargin<3
  disp('skip set to 1')
  skip=1;
end
for i=start:skip:final
  figure('PaperPosition',[0.25 2.5 34.0 12.0],'color','w','visible','off')
  subplot(1,2,1)
    SPH_plot(i)
  subplot(1,2,2)
    vortex_plot(i,'line','magnetic')
  %%%%%%%%%%%%%%%print%%%%%%%%%%%%%
  fOUT=sprintf('data/SPH%03d.png',i)
  print('-dpng',fOUT)
  close all
end
