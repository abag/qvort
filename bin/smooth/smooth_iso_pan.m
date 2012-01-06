function smooth_iso_anim(start,finish,skip) 
counter=0;
Az_new=0.;
for i=start:finish
    for j=1:5
      figure('visible','off')
      smoothed_field(i,'iso');
      counter=counter+1;
      [Az, El] = view;
      El       = ceil(El);
      Az_new = Az_new+.1;
      view([Az_new, El]);
      drawnow
      box on
      xlabel('') ; ylabel('') ; zlabel('');
      set(gca,'xtick',[]) ; set(gca,'ytick',[]) ; set(gca,'ztick',[]) ;
      fOUT=sprintf('./data/smooth_vortex%03d.png',counter)
      print('-dpng',fOUT);
      fclose('all')
    end
end
