function mesh_structure_func(x,ux,uy,uz,n,fluid)
%longitudinal structure functions - only go up to half box
sf1(1:n/2)=0.;
for k=1:n/2 ; for j=1:n ; for i=1:n
  for r=1:n/2
    sf1(r)=sf1(r)+(uz(k+r,j,i)-uz(k,j,i))^3;
  end
end ; end ; end
sf1=sf1/(n*n*n/2);
figure
plot(sf1)
%figure('Name',strcat('Iso-surface-|u|, fluid:',fluid));
%  p=patch(isosurface(x,x,x,u2,2*rms));
%  isonormals(x,x,x,u2, p)
%  set(p, 'FaceColor', 'm', 'EdgeColor', 'none');
%  daspect([1 1 1]); axis tight;
%  camup([0 0 1 ]); campos([0.7686    0.1432    0.3043])
%  camlight; lighting phong
%  xlabel('z','FontSize',14) ; ylabel('y','FontSize',14) ; zlabel('x','FontSize',14)
%  set(gca,'Fontsize',14)
%  axis([min(x) max(x) min(x) max(x) min(x) max(x)])

