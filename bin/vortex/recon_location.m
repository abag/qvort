function recon_location(t1,t2)
A=load('./data/recon_location.log');
figure('Name','recon location')
for i=1:size(A,1)
  if A(i,1)>t1 && A(i,1)<t2
     plot3(A(i,2),A(i,3),A(i,4),'ko','MarkerFaceColor','k')
     hold on
  end
end
dims=load('./data/dims.log');
axis([-dims(2)/2 dims(2)/2 -dims(2)/2 dims(2)/2 -dims(2)/2 dims(2)/2]);
box on
set(gca,'xtick',[])
set(gca,'ytick',[])
set(gca,'ztick',[])
camproj('perspective')
rotate3d on
figure('Name','recon angles')
[n xout]=histnorm(A(:,5),10);
plot(xout,n,'k','LineWidth',2)
set(gca,'FontSize',16)
xlabel('theta','FontSize',16)
ylabel('PDF','FontSize',16)
