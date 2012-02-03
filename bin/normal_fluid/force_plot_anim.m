%script is just kept to be used by batch scipt mtlb_anim.sh
%uses vortex_plot.m to create snapshots of filament
function vortex_anim(start,final,skip)
figure('visible','off');
for i=start:skip:final
  force_plot(i)
  fOUT=sprintf('data/force_plot%03d.png',i)
  print('-dpng',fOUT)
end
figure('visible','on');
close all
  
