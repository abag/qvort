%plot PDFs' (using kernal density estimates) of the x/y/z velocities along the filaments
function vortex_vel(start,finish)
if nargin<2
    finish=start
end
kde_index=6; %alter this fgor smoothness of PDF's
for i=start:finish
  filename=sprintf('data/uu%04d.dat',i);
  fid=fopen(filename);
  if fid<0
    disp('var file does not exist, exiting script')
    return
  end
  %disp(sprintf('reading var %04d',i))
  time=fread(fid,1,'float64');
  number_of_particles=fread(fid,1,'int');
  if i==start
    ux=fread(fid,number_of_particles,'float64');
    uy=fread(fid,number_of_particles,'float64');
    uz=fread(fid,number_of_particles,'float64');
  else
    l=length(ux);
    ux(l+1:l+number_of_particles)=fread(fid,number_of_particles,'float64');
    uy(l+1:l+number_of_particles)=fread(fid,number_of_particles,'float64');
    uz(l+1:l+number_of_particles)=fread(fid,number_of_particles,'float64');
  end
end
markerx=1;
markery=1;
markerz=1;
vcoff=7.
for i=1:length(ux)
    if mod(100000,i)==0
      disp(sprintf('processed %f6 percent of data',i/length(ux)))
    end
    if (ux(i)<vcoff)
        if (ux(i)>-vcoff)
            dum_ux(markerx)=ux(i);
            markerx=markerx+1;
        end
    end
    if (uy(i)<vcoff)
        if (uy(i)>-vcoff)
            dum_uy(markery)=uy(i);
            markery=markery+1;
        end
    end
     if (uz(i)<vcoff)
        if (uz(i)>-vcoff)
            dum_uz(markerz)=uz(i);
            markerz=markerz+1;
        end
    end
end       
save velocity.mat dum_ux dum_uy dum_uz
%[bwx,den_ux,xmesh]=kde(ux,2^kde_index,min(ux)*1.1,max(ux)*1.1); 
%[bwy,den_uy,ymesh]=kde(uy,2^kde_index,min(uy)*1.1,max(uy)*1.1); 
%[bwz,den_uz,zmesh]=kde(uz,2^kde_index,min(uz)*1.1,max(uz)*1.1);
%[bwz,den_uu,umesh]=kde(sqrt(ux.^2+uy.^2+uz.^2),2^kde_index,min(sqrt(ux.^2+uy.^2+uz.^2))*.8,max(sqrt(ux.^2+uy.^2+uz.^2))*1.2);
[den_ux xmesh]=hist(ux,30);
[den_uy ymesh]=hist(uy,30);
[den_uz zmesh]=hist(uz,30);
[den_uu umesh]=hist(sqrt(ux.^2+uy.^2+uz.^2),30);
histfit(ux,20) 
pause
return
subplot(2,2,1) 
  plot((xmesh),(den_ux),'-b','LineWidth',2)
  xlabel('u_x','FontSize',14)
  ylabel('PDF(u_x)','FontSize',14) 
  set(gca,'FontSize',14)
subplot(2,2,2) 
  plot(ymesh,den_uy,'-r','LineWidth',2)
  xlabel('u_y','FontSize',14)
  ylabel('PDF(u_y)','FontSize',14)
  set(gca,'FontSize',14)
subplot(2,2,3) 
  plot(zmesh,den_uz,'-m','LineWidth',2)
  xlabel('u_z','FontSize',14)
  ylabel('PDF(u_z)','FontSize',14)
  set(gca,'FontSize',14)
subplot(2,2,4) 
  plot(ymesh,den_uu,'-g','LineWidth',2)
  xlabel('|u|','FontSize',14)
  ylabel('PDF(|u|)','FontSize',14)
  set(gca,'FontSize',14)
