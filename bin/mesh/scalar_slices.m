function scalar_slices(x,v,n,fluid)
if n<=32
  interpc=4;
elseif n<=64
  interpc=2;
else
  interpc=1;
end
figure('Name',strcat('zslice-fluid: ', fluid));
  zslice(1:n,1:n)=v(n/2,:,:);
  imagesc(interp(x,interpc),interp(x,interpc),interp2(zslice,interpc))
  xlabel('x','FontSize',14) ; ylabel('y','FontSize',14)
  colorbar
  set(gca,'Fontsize',14)
figure('Name',strcat('yslice-fluid: ', fluid));
  yslice(1:n,1:n)=v(:,n/2,:);
  imagesc(interp(x,interpc),interp(x,interpc),interp2(yslice,interpc))
  xlabel('x','FontSize',14) ; ylabel('z','FontSize',14)
  colorbar
  set(gca,'Fontsize',14)
figure('Name',strcat('xslice-fluid: ', fluid));
  xslice(1:n,1:n)=v(:,:,n/2);
  imagesc(interp(x,interpc),interp(x,interpc),interp2(xslice,interpc))
  xlabel('y','FontSize',14) ; ylabel('z','FontSize',14)
  colorbar
  set(gca,'Fontsize',14)
figure('Name',strcat('Iso-surface, fluid:',fluid));
  rms=sqrt(sum(sum(sum(v.^2,3)))/n^3)
  p=patch(isosurface(x,x,x,v,1.5*rms));
  isonormals(x,x,x,v, p)
  set(p, 'FaceColor', 'm', 'EdgeColor', 'none');
  daspect([1 1 1]); axis tight; 
  camup([0 0 1 ]); campos([0.7686    0.1432    0.3043])
  camlight; lighting phong
  xlabel('z','FontSize',14) ; ylabel('y','FontSize',14) ; zlabel('x','FontSize',14)
  set(gca,'Fontsize',14)
  axis([min(x) max(x) min(x) max(x) min(x) max(x)])

