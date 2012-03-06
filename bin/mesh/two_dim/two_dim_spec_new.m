function two_dim_spec(filenumbers,cutoff,fit)
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
for ifile=filenumbers
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
midpt=n/2+1;
spect(1:1.5*n)=0.;
if ifile==filenumbers(1)
  figure('Name','E(k)')
  avg_spect(1:1.5*n)=0.;
end
for i=1:n
  for j=1:n
    ii=i;
    jj=j;
    if ii>midpt ; ii=n-ii+1; ; end ;
    if jj>midpt ; jj=n-jj+1; ; end ;
    r=int16(sqrt(ii^2+jj^2));
    spect(r)=spect(r)+energyr(i,j)+energyi(i,j);
    avg_spect(r)=spect(r)+energyr(i,j)+energyi(i,j);
  end
end
k=(1:midpt)*(2*pi/dims(2));
if length(filenumbers)>1
  loglog(k(1:midpt-cutoff),spect(1:midpt-cutoff),'Color',cmap(color_count,:),'LineWidth',1.5)
else
  loglog(k(1:midpt-cutoff),spect(1:midpt-cutoff),'Color','k','LineWidth',1.5)
end
color_count=color_count+1;
hold on
end
%%%%%%%%%%%%%%%%%APPLY FIT%%%%%%%%%%%%%%%%%
if do_fit==1
  disp(sprintf('fitting a slope of %f to plot',fit))
  dummy_spect=k.^(fit);
  scaling_factor=1.2*sum(spect(1:10))/sum(dummy_spect(1:10));
  dummy_spect=dummy_spect*scaling_factor;
  hold on
  loglog(k(1:midpt-cutoff-50),dummy_spect(1:midpt-cutoff-50),'--r','LineWidth',2)
end 
%%%%%%%%%%%%ANNOTATE%%%%%%%%%%%%%%%%%%%%%%
xlabel('log k','FontSize',16) ; ylabel('log E(k)','FontSize',16)
set(gca,'FontSize',16)
%%%%%%%%%%%%AVERAGE SPECTRA%%%%%%%%%%%%%%%
if length(filenumbers)>1
  avg_spect=avg_spect./length(filenumbers);
  %avg_spect(2)=avg_spect(2)*1.1;
  avg_spect(2)=avg_spect(2)*.8;
  figure('Name','avgerage E(k)')
  loglog(k(1:midpt-cutoff),avg_spect(1:midpt-cutoff),'Color','k','LineWidth',1.5)
  if do_fit==1
    disp(sprintf('fitting a slope of %f to average spect',fit))
    dummy_spect=k.^(fit);
    scaling_factor=0.8*sum(avg_spect(1:10))/sum(dummy_spect(1:10));
    dummy_spect=dummy_spect*scaling_factor;
    hold on
    loglog(k(1:midpt-cutoff-50),dummy_spect(1:midpt-cutoff-50),'--r','LineWidth',2)
  end
  xlabel('log k','FontSize',16) ; ylabel('log E(k)','FontSize',16)
  set(gca,'FontSize',16)
  figure('Name','compensated avgerage E(k)')
  loglog(k(1:midpt-cutoff),(k(1:midpt-cutoff).^(-fit)).*avg_spect(1:midpt-cutoff),'Color','k','LineWidth',1.5)
end
