function mesh_get_phi(ux,uy,uz,N,plotme)
ux2d=squeeze(ux(N/2,:,:));
uy2d=squeeze(uy(N/2,:,:));
phi=0.*ux2d;
phix=0.*ux2d;
for i=1:N
    phix(:,i)=sum(ux2d(:,1:i),2);
end

for j=1:N
    phi(j,:)=sum(uy2d(1:j,1),1)+phix(j,:);
end
if plotme==1
  pcolor(phi) ; shading interp
  set(gca,'FontSize',16)
end
save phi.mat phi 