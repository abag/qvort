function mesh_slices(x,ux,uy,uz,n)
u2=(ux.^2+uy.^2+uz.^2);
rms=sqrt(sum(sum(sum(u2)))/(n^3));
u2=sqrt(u2);
if n<=32
  interpc=4;
elseif n<=64
  interpc=2;
else
  interpc=1;
end
figure('Name','zslice');
  zslice(1:n,1:n)=u2(n/2,:,:);
  imagesc(interp(x,interpc),interp(x,interpc),interp2(zslice,interpc))
  xlabel('x','FontSize',14) ; ylabel('y','FontSize',14)
  set(gca,'Fontsize',14)
figure('Name','yslice');
  yslice(1:n,1:n)=u2(:,n/2,:);
  imagesc(interp(x,interpc),interp(x,interpc),interp2(yslice,interpc))
  xlabel('x','FontSize',14) ; ylabel('z','FontSize',14)
  set(gca,'Fontsize',14)
figure('Name','xslice');
  xslice(1:n,1:n)=u2(:,:,n/2);
  imagesc(interp(x,interpc),interp(x,interpc),interp2(xslice,interpc))
  xlabel('y','FontSize',14) ; ylabel('z','FontSize',14)
  set(gca,'Fontsize',14)
figure('Name','Iso-surface-|u|');
  p=patch(isosurface(x,x,x,u2));
  isonormals(x,x,x,u2, p)
  set(p, 'FaceColor', 'm', 'EdgeColor', 'none');
  daspect([1 1 1]); axis tight; 
  camup([0 0 1 ]); campos([0.7686    0.1432    0.3043])
  camlight; lighting phong
  xlabel('y','FontSize',14) ; ylabel('y','FontSize',14) ; zlabel('z','FontSize',14)
  set(gca,'Fontsize',14)