function [urms nurms]=two_dim_urms(filenumber)
filename=sprintf('./data/vel_slice_2D%03d.dat',filenumber);
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
urms=sqrt(mean(mean(usupx(:,:,:).^2+usupy(:,:,:).^2+usupz(:,:,:).^2)));
nurms=sqrt(mean(mean(unormx(:,:,:).^2+unormy(:,:,:).^2+unormz(:,:,:).^2)));
