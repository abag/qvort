function particle_shadow_anim(start,final)
figure('visible','off');
for i=start:final
  particle_shadow(i)
  fOUT=sprintf('data/par_shad%03d.png',i)
  print('-dpng',fOUT)
  close all
end
