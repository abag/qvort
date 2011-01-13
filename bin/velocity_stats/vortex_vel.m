load velocity
kde_index=4
[bwx,den_ux,xmesh]=kde(ux,2^kde_index,min(ux)*1.1,max(ux)*1.1); 
[bwy,den_uy,ymesh]=kde(uy,2^kde_index,min(uy)*1.1,max(uy)*1.1); 
[bwz,den_uz,zmesh]=kde(uz,2^kde_index,min(uz)*1.1,max(uz)*1.1); 
figure('Name','ux')
  mu=mean(ux) ; sigma=std(ux) ;
  dumx=min(ux):0.01:max(ux) ;
  dumy=normpdf(dumx,mu,sigma);
  plot((xmesh),log(den_ux),'-b','LineWidth',2)
  hold on
  plot((dumx),log(dumy),'-k','LineWidth',2)
  hold off
  %axis([min(xmesh) max(xmesh) min(log(den_ux)) max(log(den_ux))])
  xlabel('u_x','FontSize',14)
  ylabel('PDF(u_x)','FontSize',14) 
  set(gca,'FontSize',14)
figure('Name','uy')
  plot(ymesh,den_uy,'-r','LineWidth',2)
  xlabel('u_y','FontSize',14)
  ylabel('PDF(u_y)','FontSize',14)
  set(gca,'FontSize',14)
figure('Name','uz')
  plot(zmesh,den_uz,'-m','LineWidth',2)
  xlabel('u_z','FontSize',14)
  ylabel('PDF(u_z)','FontSize',14)
  set(gca,'FontSize',14)
