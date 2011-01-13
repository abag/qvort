function two_dim_slice(start,finish,fitnorm)
for i=start:finish
  ff=i;
  filename=sprintf('data/vel_slice_2D%03d.dat',i);
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
  usupx=squeeze(B(3,:,:));
  usupy=squeeze(B(4,:,:));
  usupz=squeeze(B(5,:,:));
  n=s;
  fux=fftn(usupx)/(n^2);
  fuy=fftn(usupy)/(n^2);
  fuz=fftn(usupz)/(n^2);
  energyr=real(fux).^2+real(fuy).^2+real(fuz).^2;
  energyi=imag(fux).^2+imag(fuy).^2+imag(fuz).^2;
  midpt=n/2+1;
  spect(1:1.5*n)=0.;
  for i=1:n
    for j=1:n
      ii=i;
      jj=j;
      if ii>midpt ; ii=n-ii+1; ; end ;
      if jj>midpt ; jj=n-jj+1; ; end ;
      r=int16(sqrt(ii^2+jj^2));
      spect(r)=spect(r)+energyr(i,j)+energyi(i,j);
    end
  end
  k=1:midpt;
  k2=(floor(midpt/2):midpt);
  figure('visible','off')
  subplot(2,2,1)
  loglog(k,spect(1:midpt),'LineWidth',2)
  xlabel('log k','FontSize',14) ; ylabel('log E(k)','FontSize',14)
  axis tight
  set(gca,'FontSize',14)
  dummy_spect=k.^(-5/3);
  scaling_factor=sum(spect(1:10))/sum(dummy_spect(1:10));
  dummy_spect=dummy_spect*scaling_factor;
  hold on
  loglog(k,dummy_spect,'--k','LineWidth',2)
  subplot(2,2,2)
  comp_spec=(k.^(5/3));
  loglog(k,comp_spec.*spect(1:midpt),'LineWidth',2)
  xlabel('log k','FontSize',14) ; ylabel('log E(k)*k^{5/3}','FontSize',14)
  axis tight
  set(gca,'FontSize',14)
  subplot(2,2,3)
  comp_spec=(k.^1);
  loglog(k,comp_spec.*spect(1:midpt),'LineWidth',2)
  xlabel('log k','FontSize',14) ; ylabel('log E(k)*k','FontSize',14)
  axis tight
  set(gca,'FontSize',14)
  fOUT=sprintf('sup_spec_2D%03d.png',ff);
  print('-dpng',fOUT)
  close all
  fclose(fid) ;
  clear A B ;
end
