function vortex_anim(start,final,skip)
if nargin<1
  disp('I need start and finish filenumbers')
  return
elseif nargin<2
  disp('Assuming number given is final, start set to 1')
  final=start;
  start=1;
elseif nargin<3
  disp('skip set to 1')
  skip=1;
end
for i=start:skip:final
  vortex_plot(i)
  fOUT=sprintf('data/var%03d.png',i)
  print('-dpng',fOUT)
end
  
