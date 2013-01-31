function two_dim_spec(filenumbers,cutoff,fit)
global k spect_flux
if nargin<2
    do_fit=0;
    cutoff=0;
elseif nargin<3
    do_fit=0;
else
    do_fit=1;
end
cmap=colormap(jet(length(filenumbers)));
color_count=1;
counter=0;
for ifile=filenumbers
    counter=counter+1
    filename=sprintf('data/vel_slice_2D%04d.dat',ifile);
    fid=fopen(filename);
    if fid<0
        disp('2D slice file does not exist, exiting script')
        return
    end
    dims=load('./data/dims.log');
    load data/dims.log;
    msize=dims(8);
    if (msize==0)
        disp('2D mesh size is zero exiting script')
        return
    end
    x=fread(fid,msize,'float64');
    unormx=fread(fid,msize^2,'float64');
    unormy=fread(fid,msize^2,'float64');
    unormz=fread(fid,msize^2,'float64');
    usupx=fread(fid,msize^2,'float64');
    usupy=fread(fid,msize^2,'float64');
    usupz=fread(fid,msize^2,'float64');
    unormx=reshape(unormx,msize,msize);
    unormy=reshape(unormy,msize,msize);
    unormz=reshape(unormz,msize,msize);
    usupx=reshape(usupx,msize,msize);
    usupy=reshape(usupy,msize,msize);
    usupz=reshape(usupz,msize,msize);
    xx=x;
    yy=x;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%SPECTRA%%%%%%%%%%%%%%%%%%%
    n=msize;
    fux=fftn(usupx)/(n^2);
    fuy=fftn(usupy)/(n^2);
    fuz=fftn(usupz)/(n^2);
    energyr=real(fux).^2+real(fuy).^2+real(fuz).^2;
    energyi=imag(fux).^2+imag(fuy).^2+imag(fuz).^2;
    %normalise-----------------------
    %disp('normalising')
    %energyr=energyr/sum(sum(energyr));
    %energyi=energyi/sum(sum(energyi));
    %normalise-----------------------
    midpt=n/2+1;
    spect(1:1.5*n)=0.;
    if ifile==filenumbers(1)
      spect_flux(1:length(filenumbers),1:1.5*n)=0.;
    end
    for i=1:n
        for j=1:n
            ii=i;
            jj=j;
            if ii>midpt ; ii=n-ii+1; ; end ;
            if jj>midpt ; jj=n-jj+1; ; end ;
            r=int16(sqrt(ii^2+jj^2));
            spect_flux(counter,r)=spect_flux(counter,r)+energyr(i,j)+energyi(i,j);
        end
    end
    k=(1:midpt)*(2*pi/dims(2));
end
%%%%%%%%%%%%%%%%%PLOT%%%%%%%%%%%%%%%%%
%loglog(k(1:midpt),spect_flux(:,1:midpt)')
%now compute the flux
dt=100*1E-3;
for i=1:size(spect_flux,2)
  Pdot(:,i)=gradient(spect_flux(:,i),dt);
end
plot(Pdot')
for i=1:size(Pdot,1)
intP(i,1)=0.;
i
for j=2:midpt
  intP(i,j)=trapz(Pdot(i,1:j)',k(1:j));
end
end
mean_intP=mean(intP,1);
semilogx(k,mean_intP)
hold on
semilogx(k,zeros(1,midpt))
