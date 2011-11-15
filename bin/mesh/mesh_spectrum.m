function mesh_spectrum(ux,uy,uz,n,fluid,fit)
dims=load('./data/dims.log');
fux=fftn(ux)/(n^3);
fuy=fftn(uy)/(n^3);
fuz=fftn(uz)/(n^3);
energyr=real(fux).^2+real(fuy).^2+real(fuz).^2;
energyi=imag(fux).^2+imag(fuy).^2+imag(fuz).^2;
midpt=n/2+1;
spect(1:1.5*n)=0.;
for i=1:n
    for j=1:n
        for k=1:n
            ii=i;
            jj=j;
            kk=k;
            if ii>midpt ; ii=n-ii+1; ; end ;
            if jj>midpt ; jj=n-jj+1; ; end ; 
            if kk>midpt ; kk=n-kk+1; ; end ;
            r=round(sqrt(ii^2+jj^2+kk^2));
            spect(r)=spect(r)+energyr(i,j,k)+energyi(i,j,k);
        end
    end
end
figure('Name',strcat('Energy Spectrum, fluid:',fluid)) 
k=(1:midpt)*(2*pi/dims(2));
loglog(k,spect(1:midpt),'LineWidth',2)
if fit==1
  dummy_spect=k.^(-5/3);
  scaling_factor=spect(3)/dummy_spect(3);
  dummy_spect=dummy_spect*scaling_factor;
  hold on
  loglog(k,dummy_spect,'--k','LineWidth',2)
end
xlabel('log k','FontSize',14) ; ylabel('log E(k)','FontSize',14)
axis tight
set(gca,'FontSize',14)
%save spect.mat k spect dummy_spect
