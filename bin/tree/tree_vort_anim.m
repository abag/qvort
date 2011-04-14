function tree_anim(start,final,skip)
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
figure('visible','off');
for i=start:skip:final
  subplot(1,2,1)
    octree(i)
  subplot(1,2,2)
    vortex_plot(i)
  set(gcf,'PaperUnits','centimeters')
  xSize = 30; ySize = 12; 
  xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
  set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
  fOUT=sprintf('data/octree%04d.png',i)
  print('-dpng',fOUT)
end
figure('visible','on');
close all
