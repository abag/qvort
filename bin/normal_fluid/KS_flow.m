function compressible
close all
filename='data/KS_mesh.dat';
load data/nfm_dims.log;
msize=nfm_dims(1)
if (msize==0) 
  disp('mesh size is zero exiting script')
  return
end
fid=fopen(filename);
if fid<0
  disp('mesh file does not exist, exiting script')
  return
end
x=fread(fid,msize,'float64');
ux=fread(fid,msize^3,'float64');
uy=fread(fid,msize^3,'float64');
uz=fread(fid,msize^3,'float64');
fclose(fid)
ux=reshape(ux,msize,msize,msize);
uy=reshape(uy,msize,msize,msize);
uz=reshape(uz,msize,msize,msize);
unorm=max(sqrt(ux.^2+uy.^2+uz.^2));
%plot slices of field+isosurface
mesh_slices(x,ux,uy,uz,msize,'forcing: vel field')
A=load('data/KSwavenumbers.log');
figure('Name','KS k unit vectors')
  plot3(A(:,2),A(:,3),A(:,4),'o','MarkerFaceColor','k','MarkerEdgeColor','k')
    set(gca,'FontSize',14)
    xlabel('k_x','FontSize',14)
    ylabel('k_y','FontSize',14)
    zlabel('k_z','FontSize',14)
figure('Name','KS k vectors')
  plot(A(:,1),'o','MarkerFaceColor','k','MarkerEdgeColor','k')
  set(gca,'FontSize',14)
    xlabel('N','FontSize',14)
    ylabel('k_N','FontSize',14)
  