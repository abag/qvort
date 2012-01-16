function mesh_plot(start,finish,skip)
load data/sm_dims.log;
msize=sm_dims(1);
counter=1;
for i=start:skip:finish
filename=sprintf('data/smoothed_field%03d.dat',i);
fid=fopen(filename);
if fid<0
  disp('file does not exist, exiting script')
  return
end
t=fread(fid,1,'float64');
x=fread(fid,msize,'float64');
wx=fread(fid,msize^3,'float64');
wy=fread(fid,msize^3,'float64');
wz=fread(fid,msize^3,'float64');
fclose(fid);
rms=sqrt(mean(wx.^2+wy.^2+wz.^2));
tt(counter)=t;
w_rms(counter)=rms;
counter=counter+1;
end
plot(tt,w_rms,'-k','LineWidth',2)
xlabel('t','FontSize',16)
ylabel('\Omega_{rms}','FontSize',16)
set(gca,'FontSize',16)
