function povray_anim(start,final,skip)
if nargin==2
    skip=1
end
for i=start:skip:final
  vortex_load(i)
  povray_out
  system(strcat('/opt/local/bin/povray +I./povray/povray_plot.pov +O',sprintf('./pov_snap%03d.jpg',i),' +Fj +W3200 +H2400 +V -D +X'))
end
delete('./test.log')
delete('./test.loge')
