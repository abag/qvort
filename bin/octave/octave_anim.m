function octave_anim(start,final,skip)
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
  octave_plot(i)
  fOUT=sprintf('data/var%04d.png',i)
  print('-dpng',fOUT)
end
figure('visible','on');
close all
  
