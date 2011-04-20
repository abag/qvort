%script is just kept to be used by batch scipt mtlb_anim.sh
%uses vortex_plot.m to create snapshots of filament
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
figure('visible','off');
for i=start:skip:final
  vortex_plot(i,'line','magnetic')
  fOUT=sprintf('data/var%04d.eps',i)
  print('-depsc',fOUT)
end
figure('visible','on');
close all
  
