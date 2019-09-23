%script is just kept to be used by batch scipt mtlb_anim.sh
%uses vortex_plot.m to create snapshots of filament
function vortex_anim(start,final,skip,varargin)
eps=0;
figure('visible','off');
for i=start:skip:final
  %vortex_plot(i,'line','magnetic')
  vortex_plot(i,varargin{:})
  %vortex_plot_sep_loops(i)
  %view(0.2*i,26)
  if eps==1 
    fOUT=sprintf('data/var%04d.eps',i)
    print('-depsc',fOUT)
  else
    fOUT=sprintf('data/var%04d.jpeg',i)
    print('-djpeg',fOUT)
  end 
  %close all open files
  fclose('all');
end
figure('visible','on');
close all
  
