function one_dim_lattice(filenumber,fit)
if nargin<2
  do_fit=0;
else
  do_fit=1;
end
plot_slice=0;
plot_lines=0;
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
  if plot_slice==1
    plot(meshus_sq(i,:))
    pause
  end
end
%scale velocity into a colormap
if plot_lines==1
  store_caxis=([min(min(log(meshus_sq))) max(max(log(meshus_sq)))]);
  meshus_sq=log(meshus_sq)-min(min(log(meshus_sq)));
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
end
%SPECTRA
save out.mat meshus
for i=1:nlines
  fux=fftn(squeeze(meshus(i,:,1)))/nmesh;
  fuy=fftn(squeeze(meshus(i,:,2)))/nmesh;
  fuz=fftn(squeeze(meshus(i,:,3)))/nmesh;
  energyr=real(fux).^2+real(fuy).^2+real(fuz).^2;
  energyi=imag(fux).^2+imag(fuy).^2+imag(fuz).^2;
  midpt=nmesh/2+1;
  spect(1:1.5*nmesh)=0.;
  for j=1:nmesh
    jj=j;
    if jj>midpt ; jj=nmesh-jj+1; ; end ;
    spec(jj)=spect(jj)+energyr(j)+energyi(j);
  end
  k=1:midpt-1;
  loglog(k,spect(1:midpt-1),'LineWidth',2)
  set(gca,'FontSize',14)
  ylabel('log E(k)','FontSize',14)
  xlabel('log k','FontSize',14)
  if do_fit==1
    disp(sprintf('fitting a slope of %f to plot',fit))
    dummy_spect=k.^(fit);
    scaling_factor=sum(spect(1:10))/sum(dummy_spect(1:10));
    dummy_spect=dummy_spect*scaling_factor;
    hold on
    loglog(k,dummy_spect,'--k','LineWidth',2)
    hold off
  end
  pause
end

