function mesh_slices(x,ux,uy,uz,n,fluid)
u2=(ux.^2+uy.^2+uz.^2);
u2=sqrt(u2);
if n<=32
  interpc=4;
elseif n<=64
  interpc=2;
else
  interpc=1;
end
figure('Name',strcat('zslice-fluid: ', fluid));
  zslice(1:n,1:n)=u2(n/2,:,:);
  imagesc(interp(x,interpc),interp(x,interpc),interp2(zslice,interpc))
  set(gca,'YDir','normal')
  hold on
  quiver(x,x,squeeze(ux(n/2,:,:)),squeeze(uy(n/2,:,:)),'k')
  xlabel('x','FontSize',14) ; ylabel('y','FontSize',14)
  colorbar
  set(gca,'Fontsize',14)
figure('Name',strcat('yslice-fluid: ', fluid));
  yslice(1:n,1:n)=u2(:,n/2,:);
  imagesc(interp(x,interpc),interp(x,interpc),interp2(yslice,interpc))
  set(gca,'YDir','normal')
  xlabel('x','FontSize',14) ; ylabel('z','FontSize',14)
  colorbar
  set(gca,'Fontsize',14)
figure('Name',strcat('xslice-fluid: ', fluid));
  xslice(1:n,1:n)=u2(:,:,n/2);
  imagesc(interp(x,interpc),interp(x,interpc),interp2(xslice,interpc))
  set(gca,'YDir','normal')
  xlabel('y','FontSize',14) ; ylabel('z','FontSize',14)
  colorbar
  set(gca,'Fontsize',14)

