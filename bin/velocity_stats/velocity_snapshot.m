function velocity_snapshot(filenumber)
  filename=sprintf('data/uu%04d.dat',filenumber);
  fid=fopen(filename);
  if fid<0
    disp('var file does not exist, exiting script')
    return
  end
  %disp(sprintf('reading var %04d',i))
  time=fread(fid,1,'float64');
  number_of_particles=fread(fid,1,'int');
  ux=fread(fid,number_of_particles,'float64');
  uy=fread(fid,number_of_particles,'float64');
  uz=fread(fid,number_of_particles,'float64');
  subplot(1,3,1)
    [f,xi]=ksdensity(ux);
    plot(xi,f,'r-','LineWidth',2)
    xlabel('u_x','FontSize',14)
    ylabel('PDF','FontSize',14)
    set(gca,'FontSize',14)
    title(strcat('mean= ',num2str(mean(ux))),'FontSize',14)
  subplot(1,3,2)
    [f,xi]=ksdensity(uy);
    plot(xi,f,'g-','LineWidth',2)
    xlabel('u_y','FontSize',14)
    ylabel('PDF','FontSize',14)
    title(strcat('mean= ',num2str(mean(uy))),'FontSize',14)
  subplot(1,3,3)
    [f,xi]=ksdensity(uz);
    plot(xi,f,'b-','LineWidth',2)
    xlabel('u_z','FontSize',14)
    ylabel('PDF','FontSize',14)
    title(strcat('mean= ',num2str(mean(uz))),'FontSize',14)
end
