%not ready yet
function vortex_fft(filenumber)
filename=sprintf('data/var%03d.log',filenumber);
%set options based on varargin
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%get the dimensions information from dims.log
dims=load('./data/dims.log');
if dims(4)==1
  fid=fopen(filename);
  if fid<0
      disp('var file does not exist, exiting script')
      return
  end
  time=fread(fid,1,'float64');
  number_of_particles=fread(fid,1,'int');
  x=fread(fid,number_of_particles,'float64');
  y=fread(fid,number_of_particles,'float64');
  z=fread(fid,number_of_particles,'float64');
  f=fread(fid,number_of_particles,'int');
  u=fread(fid,number_of_particles,'float64');
else 
  fid=fopen(filename);
  if fid<0
      disp('var file does not exist, exiting script')
      return
  end
  %read the time
  tline=fgetl(fid);
  dummy=textscan(tline, '%f');
  time=dummy{:};
  %how many particles
  tline=fgetl(fid);
  dummy=textscan(tline, '%d');  
  number_of_particles=dummy{:};
  %get the particles%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  for j=1:number_of_particles
    tline=fgetl(fid);
    dummy=textscan(tline, '%f');
    dummy_vect=dummy{:};
    x(j)=dummy_vect(1);
    y(j)=dummy_vect(2);
    z(j)=dummy_vect(3);
    f(j)=dummy_vect(4);
    u(j)=dummy_vect(5);
  end
  f=uint16(f);
end
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
            r=int16(sqrt(ii^2+jj^2+kk^2));
            spect(r)=spect(r)+energyr(i,j,k)+energyi(i,j,k);
        end
    end
end
figure('Name',strcat('Energy Spectrum, fluid:',fluid)) 
k=1:midpt;
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
