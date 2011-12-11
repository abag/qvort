function two_dim_spec(filenumber,order)
if nargin<2
  order=3;
end
filename=sprintf('data/vel_slice_2D%04d.dat',filenumber);
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
fclose(fid)
%%%%%%%%%%%%%%%%%%%%STRUCTURE FUNCTION%%%%%%%%%%%%%%%%%%%
n=s;
sf(1:n/2)=0.
for i=1:n/2 ; for j=1:n
  for r=1:n/2
    sf(r)=sf(r)+(usupy(i+r,j)-usupy(i,j))^pow;
  end
end ; end
for i=1:n/2 ; for j=1:n
  for r=1:n/2
    sf(r)=sf(r)+(usupy(j,i+r)-usupy(j,i))^pow;
  end
end ; end
%%%%%%%%%%%%ANNOTATE%%%%%%%%%%%%%%%%%%%%%%
plot([1:n/2],sf)
xlabel('log k','FontSize',16) ; ylabel('log E(k)','FontSize',16)
set(gca,'FontSize',16)
