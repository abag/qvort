function mesh_double_spectrum(ux,uy,uz,n,fluid,fit)
dims=load('./data/dims.log');
uu=sqrt(ux.^2+uy.^2+uz.^2);
vcoff=5000. ;
index = find(uu > vcoff);
ux(index)=0. ; uy(index)=0. ;uz(index)=0. ;
clear index
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
            ii=i-1;
            jj=j-1;
            kk=k-1;
            if ii>=midpt ; ii=n-ii; ; end ;
            if jj>=midpt ; jj=n-jj; ; end ; 
            if kk>=midpt ; kk=n-kk; ; end ;
            r=round(sqrt(ii^2+jj^2+kk^2));
            if r==0 ; continue ; end
            spect(r)=spect(r)+energyr(i,j,k)+energyi(i,j,k);
        end
    end
end
figure('Name','Normal (blue) and super (black) fluid energy spectra') 
k=(1:midpt)*(2*pi/dims(2));
loglog(k,spect(1:midpt),'LineWidth',2)
dummy_spect=k.^(-5/3);
scaling_factor=spect(3)/dummy_spect(3);
dummy_spect=dummy_spect*scaling_factor;
hold on
loglog(k,dummy_spect,'--r','LineWidth',2)
xlabel('log k','FontSize',14) ; ylabel('log E(k)','FontSize',14)
axis tight
set(gca,'FontSize',14)
