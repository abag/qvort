function full_vd_within_smooth(start,final,skip)
smooth_factor=2;
for i=start:skip:final
  [t(i) ld(i)]=vortex_density_within_smooth(i,smooth_factor*i); 
end
plot(t,ld,'k-','LineWidth',2)
set(gca,'FontSize',16)
xlabel('t','FontSize',16)
ylabel('L_{smooth}','FontSize',16)
save density_within_smooth.mat t ld 
  
