function two_dim_spec_new_fitting(filenumbers)
do_fit=1;
cutoff=0;
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
intervortex=mean(1./sqrt(A(floor(0.5*length(A)):length(A),6)/dims(2)^3));
%%%%%%%%%%%%%CHECKING FOR ABNORMAL VELOCITY%%%%%%%%%%%%%%%%
uu=sqrt(usupx.^2+usupy.^2+usupz.^2);
vcoff=50. ;
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
    ii=i;
    jj=j;
    if ii>midpt ; ii=n-ii+1; ; end ;
    if jj>midpt ; jj=n-jj+1; ; end ;
    r=int16(sqrt(ii^2+jj^2));
    spect(r)=spect(r)+energyr(i,j)+energyi(i,j);
    avg_spect(r)=avg_spect(r)+energyr(i,j)+energyi(i,j);
  end
end
k=(1:midpt)*(2*pi/dims(2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if length(filenumbers)>1
%   loglog(k(1:midpt-cutoff),spect(1:midpt-cutoff),'Color',cmap(color_count,:),'LineWidth',1.5)
% else
%   loglog(k(1:midpt-cutoff),spect(1:midpt-cutoff),'Color','k','LineWidth',1.5)
% end
% color_count=color_count+1;
% hold on
 end
%%%%%%%%%%%%%%%%%APPLY FIT%%%%%%%%%%%%%%%%%
% if do_fit==1
%   disp(sprintf('fitting a slope of %f to plot',fit))
%   dummy_spect=k.^(fit);
%   scaling_factor=1.2*sum(spect(1:10))/sum(dummy_spect(1:10));
%   dummy_spect=dummy_spect*scaling_factor;
%   hold on
%   loglog(k(midpt-cutoff-80:midpt-cutoff),dummy_spect(midpt-cutoff-80:midpt-cutoff),'--r','LineWidth',2)
% end 
%%%%%%%%%%%%ANNOTATE%%%%%%%%%%%%%%%%%%%%%%
% xlabel('log k','FontSize',16) ; ylabel('log E(k)','FontSize',16)
% set(gca,'FontSize',16)
%%%%%%%%%%%%AVERAGE SPECTRA%%%%%%%%%%%%%%%
if length(filenumbers)>=1
  avg_spect=avg_spect./length(filenumbers);
  %avg_spect(2)=avg_spect(2)*1.1;
  avg_spect(2)=avg_spect(2)*.8;
  k_box = round(dims(2)/dims(2));
  k_intervortex = round(dims(2)/intervortex);
  k_resolution = round(dims(2)/dims(1));
  Ec = trapz(k(k_box:k_intervortex),avg_spect(k_box:k_intervortex));
  Eq = trapz(k(k_intervortex:k_resolution),avg_spect(k_intervortex:k_resolution));
  disp(sprintf('Ec/Eq is: '))
  Ec/Eq
  figure('Name','avgerage E(k)')
  loglog(k(1:k_resolution),avg_spect(1:k_resolution),'Color','k','LineWidth',1.5)
  if do_fit==1
    disp(sprintf('I have estimated the intervortex spacing from the final half of the data, l=%f',intervortex))
    disp(sprintf('creating a fit from k_intervortex to k_resolution'))
    cfit = polyfit(log(k(k_intervortex:k_resolution)),log(avg_spect(k_intervortex:k_resolution)),1);
    dummy_spect=k.^(cfit(1));
    disp(sprintf('the fit values are: '))
    cfit
    scaling_factor=0.8*sum(avg_spect(1:10))/sum(dummy_spect(1:10));
    dummy_spect=dummy_spect*scaling_factor;
    hold on
    loglog(k(k_intervortex:k_resolution),dummy_spect(k_intervortex:k_resolution),'--r','LineWidth',2)
  end
  xlabel('log k','FontSize',16) ; ylabel('log E(k)','FontSize',16)
  set(gca,'FontSize',16)
  figure('Name','compensated avgerage E(k)')
  loglog(k(k_intervortex:k_resolution),(k(k_intervortex:k_resolution).^(-cfit(1))).*avg_spect(k_intervortex:k_resolution),'Color','k','LineWidth',1.5)
end
