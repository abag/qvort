function two_dim_spec_new(filenumbers,fit)
if nargin<2
  do_fit=0;
else
  do_fit=1;
end
cmap=colormap(jet(length(filenumbers)));
color_count=1;
for ifile=filenumbers
filename=sprintf('./data/vel_slice_2D%04d.dat',ifile);
fid=fopen(filename);
if fid<0
  disp('2D slice file does not exist, exiting script')
  return
end
dims=load('./data/dims.log');
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
A=load('./data/ts.log');
if length(filenumbers)>1
    %intervortex=mean(1./sqrt(A(floor(0.5*length(A)):length(A),6)/dims(2)^3));
end
%%%%%%%%%%%%%CHECKING FOR ABNORMAL VELOCITY%%%%%%%%%%%%%%%%
uu=sqrt(usupx.^2+usupy.^2+usupz.^2);
vcoff=10. ;
index = find(uu > vcoff);
usupx(index)=0. ; usupy(index)=0. ;usupz(index)=0. ;
uu=sqrt(usupx.^2+usupy.^2+usupz.^2);
clear index
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
        ii=i-1;
        jj=j-1;
        if ii>=midpt ; ii=n-ii;  end ;
        if jj>=midpt ; jj=n-jj;  end ;
        r=round(sqrt(ii^2+jj^2));
        if r==0 ; continue ; end
        spect(r)=spect(r)+energyr(i,j)+energyi(i,j);
    end
end
avg_spect=avg_spect+spect;
k=(1:midpt)*(2*pi/dims(2));
k_resolution = round(dims(2)/dims(1));
if length(filenumbers)>1
  loglog(k(1:k_resolution),spect(1:k_resolution),'Color',cmap(color_count,:),'LineWidth',1.5)
else
  loglog(k(1:k_resolution),spect(1:k_resolution),'Color','k','LineWidth',1.5)
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
  loglog(k(1:k_resolution-50),dummy_spect(1:k_resolution-50),'--r','LineWidth',2)
end 
%%%%%%%%%%%%ANNOTATE%%%%%%%%%%%%%%%%%%%%%%
xlabel('log k','FontSize',16) ; ylabel('log E(k)','FontSize',16)
set(gca,'FontSize',16)
%%%%%%%%%%%%AVERAGE SPECTRA%%%%%%%%%%%%%%%
if length(filenumbers)>1
  avg_spect=avg_spect./length(filenumbers);
  %avg_spect(2)=avg_spect(2)*1.1;
  figure('Name','avgerage E(k)')
  loglog(k(1:k_resolution),avg_spect(1:k_resolution),'Color','k','LineWidth',1.5)
  if do_fit==1
    disp(sprintf('fitting a slope of %f to average spect',fit))
    dummy_spect=k.^(fit);
    scaling_factor=0.8*sum(avg_spect(1:10))/sum(dummy_spect(1:10));
    dummy_spect=dummy_spect*scaling_factor;
    hold on
    k_intervortex = round(dims(2)/intervortex);
    loglog(k(1:k_intervortex),dummy_spect(1:k_intervortex),'--r','LineWidth',2)
  end
  xlabel('log k','FontSize',16) ; ylabel('log E(k)','FontSize',16)
  set(gca,'FontSize',16)
  figure('Name','compensated avgerage E(k)')
  loglog(k(1:k_resolution),(k(1:k_resolution).^(-fit)).*avg_spect(1:k_resolution),'Color','k','LineWidth',1.5)
  compensated_spect=k(1:k_resolution).^(-fit).*avg_spect(1:k_resolution);
end
save spec.mat k avg_spect k_resolution
if length(filenumbers)>1
    save compensated_spec.mat k k_resolution compensated_spect
end
