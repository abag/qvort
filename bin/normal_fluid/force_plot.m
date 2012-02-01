function force_plot(filenumber)
close all
filename=sprintf('./data/forcing_mesh%03d.dat',filenumber);
msize=64;
fid=fopen(filename);
if fid<0
  disp('mesh file does not exist, exiting script')
  return
end
x=fread(fid,msize,'float64');
ux=fread(fid,msize^2,'float64');
uy=fread(fid,msize^2,'float64');
uz=fread(fid,msize^2,'float64');
fclose(fid)
ux=reshape(ux,msize,msize);
uy=reshape(uy,msize,msize);
uz=reshape(uz,msize,msize);
unorm=(sqrt(ux.^2+uy.^2+uz.^2));
pcolor(x,x,unorm) ; shading interp
set(gca,'FontSize',16)
xlabel('x','FontSize',16)
ylabel('y','FontSize',16)
colormap(hot)
hold on
quiver(x,x,ux,uy,'k')
axis equal
axis tight

