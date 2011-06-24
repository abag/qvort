function mesh_iso(x,ux,uy,uz,n,fluid)
u2=(ux.^2+uy.^2+uz.^2);
rms=sqrt(sum(sum(sum(u2)))/(n^3))
u2=sqrt(u2);
figure('Name',strcat('Iso-surface-|u|, fluid:',fluid));
  p=patch(isosurface(x,x,x,u2,rms));
  isonormals(x,x,x,u2, p)
  set(p, 'FaceColor', 'm', 'EdgeColor', 'none');
  daspect([1 1 1]); axis tight;
  camup([0 0 1 ]); campos([0.7686    0.1432    0.3043])
  camlight; lighting phong
  xlabel('z','FontSize',14) ; ylabel('y','FontSize',14) ; zlabel('x','FontSize',14)
  set(gca,'Fontsize',14)
  axis([min(x) max(x) min(x) max(x) min(x) max(x)])

