function smooth_iso_anim(start,finish,skip) 
counter=0;
campan_amount=0;
for i=start:finish
    for j=1:10
      smoothed_field(i,'iso')
      counter=counter+1;
      campan(campan_amount,0.);
      campan_amount=campan_amount+18;
      if campan_amount==360
        campan_amount=0 ; %reset view
      end 
    end
    fOUT=sprintf('./data/smooth_vortex%03d.png',counter);
    print('-dpng',fOUT);
end
