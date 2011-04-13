function mesh_structure_func(x,ux,uy,uz,n,fluid)
%longitudinal structure functions - only go up to half box
sf1(1:n/2)=0.;
pow=2;
for k=1:n/2 ; for j=1:n ; for i=1:n
  for r=1:n/2
    sf1(r)=sf1(r)+abs(uz(k+r,j,i)-uz(k,j,i))^pow;
  end
end ; end ; end
for k=1:n ; for j=1:n/2 ; for i=1:n
  for r=1:n/2
    sf1(r)=sf1(r)+abs(uy(k,j+r,i)-uy(k,j,i))^pow;
  end
end ; end ; end
for k=1:n ; for j=1:n ; for i=1:n/2
  for r=1:n/2
    sf1(r)=sf1(r)+abs(ux(k,j,i+r)-ux(k,j,i))^pow;
  end
end ; end ; end
sf1=sf1/(3*n*n*n/2);
figure('Name',strcat('2nd order struct func., fluid:',fluid));
loglog((1:n/2),(sf1(1:n/2)),'LineWidth',2)
xlabel('r','FontSize',14)
ylabel('struct. func.','FontSize',14) 
set(gca,'Fontsize',14)

