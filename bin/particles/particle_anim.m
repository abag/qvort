function particle_anim(start,final,skip)
vortex=0 %if set to 1 plots vortices
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
  if vortex==1
    vortex_plot(i)
    hold on
  end
  particle_plot(i)
  fOUT=sprintf('data/par%03d.jpeg',i)
  print('-djpeg',fOUT)
  close all
end
