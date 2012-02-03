%script is just kept to be used by batch scipt mtlb_anim.sh
%uses vortex_point_slice.m to create snapshots of point slices in x, y or z planes
function vortex_point_slice_anim(start,final,skip,plane)
eps=0;
figure('visible','off');
for i=start:skip:final
  vortex_point_slice(i,plane)
  if eps==1 
    fOUT=sprintf('data/var%04d.eps',i)
    print('-depsc',fOUT)
  else
    fOUT=sprintf('data/var%04d.png',i)
    print('-dpng',fOUT)
  end 
end
figure('visible','on');
close all
  
