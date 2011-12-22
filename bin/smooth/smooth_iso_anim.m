function smooth_iso_anim(start,finish,skip) 
for i=start:finish
    smoothed_field(i,'iso')
    fOUT=sprintf('smooth_vortex%03d.png',i);
    print('-dpng',fOUT);
end
