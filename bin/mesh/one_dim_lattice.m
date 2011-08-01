function one_dim_lattice(filenumber,fit)
if nargin<2
  do_fit=0;
else
  do_fit=1;
end
%load in dims
dims=load('data/one_dim_lattice_dims.log');
nlines=dims(1);
nmesh=dims(2);
filename=sprintf('data/1D_lattice%04d.dat',filenumber);
fid=fopen(filename);
if fid<0
  disp('1D lattice file does not exist, exiting script')
  return
end
A=fread(fid,'float64');
A=reshape(A,9,nmesh,nlines);
meshx(1:nlines,1:nmesh,1:3)=0.;
meshus(1:nlines,1:nmesh,1:3)=0.;
meshus_sq(1:nlines,1:nmesh)=0.;
meshun(1:nlines,1:nmesh,1:3)=0.;
for i=1:nlines
  for j=1:nmesh
    meshx(i,j,1:3)=(A(1:3,j,i));
    meshus(i,j,1:3)=squeeze(A(4:6,j,i));
    meshus_sq(i,j)=squeeze(sqrt(A(4,j,i).^2+A(5,j,i).^2+A(6,j,i).^2));
    meshun(i,j,1:3)=squeeze(A(7:9,j,i));
  end
end
%scale velocity into a colormap
store_caxis=([min(min((meshus_sq))) max(max(meshus_sq))]);
meshus_sq=meshus_sq-min(min((meshus_sq)));
rainbow_scale=199/max(max(meshus_sq)) ;
meshus_sq=meshus_sq*rainbow_scale;
rainbowcmap=colormap(jet(200));
for i=1:nlines
  for j=1:nmesh-1
    plot3([meshx(i,j,1) meshx(i,j+1,1)], [meshx(i,j,2) meshx(i,j+1,2)], [meshx(i,j,3) meshx(i,j+1,3)],'Color',rainbowcmap(max(1,ceil(meshus_sq(i,j))),:))
    hold on
  end
end
dims2=load('./data/dims.log');
axis([-dims2(2)/2 dims2(2)/2 -dims2(2)/(2*dims2(7)) dims2(2)/(2*dims2(7)) -dims2(2)/(2*dims2(7)) dims2(2)/(2*dims2(7))]);
caxis(store_caxis)
colorbar
