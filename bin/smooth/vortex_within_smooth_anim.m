%script is just kept to be used by batch scipt mtlb_anim.sh
%uses vortex_plot.m to create snapshots of filament
function vortex_within_smooth_anim(start,final,skip,smooth_file)
eps=0;
figure('visible','off');
for i=start:skip:final
  plot_vortex_within_smooth(smooth_file,i)
  if eps==1 
    fOUT=sprintf('data/vortex_in_smooth%04d.eps',i)
    print('-depsc',fOUT)
  else
    fOUT=sprintf('data/vortex_in_smooth%04d.png',i)
    print('-dpng',fOUT)
  end 
end
figure('visible','on');
close all
  
