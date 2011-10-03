function two_dim_spec(filenumber,breakpoint)
filename=sprintf('data/vel_slice_2D%04d.dat',filenumber);
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%SPECTRA%%%%%%%%%%%%%%%%%%%
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
plot(k(1:midpt),spect(1:midpt),'k','LineWidth',2)
for i=1:length(k)-1
  if (k(i)<=breakpoint) && (k(i+1)>breakpoint)
    breakpoint_index=i;
  end
end
I1=simpsons(spect(1:breakpoint_index),k(1),k(breakpoint_index),[]) 
I2=simpsons(spect(breakpoint_index+1:400),k(breakpoint_index+1),k(400),[]) 
I1/I2

