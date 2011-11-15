function two_dim_spec(filenumbers,cutoff,fit)
if nargin<2
    do_fit=0;
    cutoff=0;
elseif nargin<3
    do_fit=0;
else
    do_fit=1;
end
for i=filenumbers
    filename=sprintf('data/vel_slice_2D%04d.dat',i);
    fid=fopen(filename);
    if fid<0
        disp('2D slice file does not exist, exiting script')
        return
    end
    dims=load('./data/dims.log');
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
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%SPECTRA%%%%%%%%%%%%%%%%%%%
    n=s;
    fux=fftn(usupx)/(n^2);
    fuy=fftn(usupy)/(n^2);
    fuz=fftn(usupz)/(n^2);
    energyr=real(fux).^2+real(fuy).^2+real(fuz).^2;
    energyi=imag(fux).^2+imag(fuy).^2+imag(fuz).^2;
    midpt=n/2+1;
    dspect(1:1.5*n)=0.;
    for i=1:n
        for j=1:n
            ii=i;
            jj=j;
            if ii>midpt ; ii=n-ii+1; ; end ;
            if jj>midpt ; jj=n-jj+1; ; end ;
            r=int16(sqrt(ii^2+jj^2));
            dspect(r)=dspect(r)+energyr(i,j)+energyi(i,j);
        end
    end
    if exist('spect')==1
       spect(:)=spect(:)+dspect(:);
    else
        spect=dspect;
    end
end
k=(1:midpt)*(2*pi/dims(2));
spect(:)=spect(:)/length(filenumbers);
loglog(k(1:midpt-cutoff),spect(1:midpt-cutoff),'k','LineWidth',2)
%%%%%%%%%%%%%%%%%APPLY FIT%%%%%%%%%%%%%%%%%
if do_fit==1
    disp(sprintf('fitting a slope of %f to plot',fit))
    dummy_spect=k.^(fit);
    scaling_factor=0.8*sum(spect(1:10))/sum(dummy_spect(1:10));
    dummy_spect=dummy_spect*scaling_factor;
    hold on
    loglog(k,dummy_spect,'--r','LineWidth',2)
end
%%%%%%%%%%%%ANNOTATE%%%%%%%%%%%%%%%%%%%%%%
xlabel('log k','FontSize',16) ; ylabel('log E(k)','FontSize',16)
set(gca,'FontSize',16)
