function two_dim_slice(start,finish)
for i=start:finish
  disp(sprintf('processing file %04d of %04d',i,finish)) 
  filename=sprintf('data/vel_slice_2D%04d.dat',i);
  fid=fopen(filename);
  if fid<0
    disp('2D slice file does not exist, exiting script')
    return
  end
  A=fread(fid,'float64');
  s=length(A); s=s/8; s=sqrt(s);
  B=reshape(A,8,s,s);
  x=squeeze(B(1,:,:));
  y=squeeze(B(2,:,:));
  xx=squeeze(x(1,1,:));
  yy=squeeze(y(1,:,1));
  usupx=squeeze(B(3,:,:));
  usupy=squeeze(B(4,:,:));
  usupz=squeeze(B(5,:,:));
  unormx=squeeze(B(6,:,:));
  unormy=squeeze(B(6,:,:));
  unormz=squeeze(B(6,:,:));
  figure('visible','off')
  uu=sqrt(usupx.^2+usupy.^2+usupz.^2);
  vcoff=1. ;
  index = find(uu > vcoff);
  uu(index) = vcoff;
  clear index
  imagesc(xx,yy,uu) ; shading interp
  set(gca,'FontSize',14)
  xlabel('x','FontSize',14)
  ylabel('y','FontSize',14)
  colorbar
  fOUT=sprintf('sup_slice_2D%03d.png',i);
  print('-dpng',fOUT)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if std(sqrt(unormx.^2+unormy.^2+unormz.^2))>0
    figure('visible','off')
    imagesc(xx,yy,sqrt(unormx.^2+unormy.^2+unormz.^2)) ; shading interp
    set(gca,'FontSize',14)
    xlabel('x','FontSize',14)
    ylabel('y','FontSize',14)
    colorbar
    fOUT=sprintf('norm_slice_2D%03d.png',i);
    print('-dpng',fOUT)
  end
end
