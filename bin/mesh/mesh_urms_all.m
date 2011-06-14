function mesh_urms_all(start,finish,skip)
counter=1;
for i=start:skip:finish
  [t urms nurms]=mesh_urms(i);
  tt(i)=t;
  u2(i)=urms;
  nu2(i)=nurms;
  counter=counter+1;
end
plot(tt,u2,'k','LineWidth',2)
set(gca,'FontSize',14)
xlabel('t','FontSize',14)
ylabel('E','FontSize',14)
