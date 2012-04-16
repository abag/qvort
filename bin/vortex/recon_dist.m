function recon_dist(filenumber)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename=sprintf('./data/recon_dist%04d.log',filenumber);
fid=fopen(filename);
if fid<0
  disp('file does not exist, exiting script')
  return
end
%read the time
tline=fgetl(fid);
dummy=textscan(tline, '%f');
time=dummy{:}
tline=fgetl(fid);
dummy=textscan(tline, '%f');
dt=dummy{:}
tline=fgetl(fid);
dummy=textscan(tline, '%f');
ndist=dummy{:}
tline=fgetl(fid);
dummy=textscan(tline, '%f');
angle=dummy{:}
for j=1:ndist
  tline=fgetl(fid);
  dummy=textscan(tline, '%f');
  dummy_vect=dummy{:};
  dist(j)=dummy_vect(1,:);
  curvi(j)=dummy_vect(2,:);
  curvj(j)=dummy_vect(3,:);
end
t=(0:ndist-1)*dt;
figure('name','natural scale')
plot(t,dist.^2,'k-','LineWidth',2)
hold on
plot(t,(2*pi*9.97E-4)*t,'r')
hold off
set(gca,'FontSize',16)
xlabel('t','FontSize',16)
ylabel('\delta^2','FontSize',16)
figure('name','log scale')
dist=dist-min(dist);
loglog(t,dist,'ko-','LineWidth',2)
hold on
loglog(t,((9.97E-4)*t).^(1/2),'r')
hold off
set(gca,'FontSize',16)
xlabel('t','FontSize',16)
ylabel('\delta','FontSize',16)
figure('name',strcat('curvature - initial angle - ',num2str(angle)))
plot(t,curvi,'-b','LineWidth',2)
hold on
plot(t,curvj,'-r','LineWidth',2)
hold off
set(gca,'FontSize',16)
xlabel('t','FontSize',16)
ylabel('curvature','FontSize',16)
t2=t(20:end);
figure;
dist2=dist(20:end);
loglog(t2,dist2)
polyfit(t2,dist2,3)



