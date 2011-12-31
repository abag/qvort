function smooth_iso_anim(start,finish,skip) 
counter=0;
for i=start:finish
    for j=1:10
      smoothed_field(i,'iso')
      counter=counter+1;
      [Az, El] = view;
      El       = ceil(El);
      Az_new = Az+18
      view([Az_new, El]);
      drawnow
    end
    fOUT=sprintf('./data/smooth_vortex%03d.png',counter);
    print('-dpng',fOUT);
end
