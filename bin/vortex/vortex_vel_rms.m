function vortex_vel_rms(start,finish)
global u u2 f
counter=0
for i=start:finish
    counter=counter+1
vortex_load(i)
vrms(counter)=mean(u2.^2)
end
end 
