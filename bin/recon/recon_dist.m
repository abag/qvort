function [a_fit b_fit c_fit R_fit angle cmax]=recon_dist(filenumber,plotme)
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
time=dummy{:};
tline=fgetl(fid);
dummy=textscan(tline, '%f');
dt=dummy{:};
tline=fgetl(fid);
dummy=textscan(tline, '%f');
ndist=dummy{:};
tline=fgetl(fid);
dummy=textscan(tline, '%f');
angle=dummy{:};
for j=1:ndist
    tline=fgetl(fid);
    dummy=textscan(tline, '%f');
    dummy_vect=dummy{:};
    dist(j)=dummy_vect(1,:);
    curvi(j)=dummy_vect(2,:);
    curvj(j)=dummy_vect(3,:);
end
%maximum curvature averaged over the two points
cmax=mean([max(curvi) max(curvj)])
t=(0:ndist-1)*dt;
if plotme>0
    figure('name','natural scale')
    plot(t,dist.^2,'k-','LineWidth',2)
    hold on
    plot(t,(2*pi*9.97E-4)*t,'r')
    hold off
    set(gca,'FontSize',16)
    xlabel('t','FontSize',16)
    ylabel('\delta^2','FontSize',16)
end
%-----------Lucy fitting----------------
t_fit=t(200:end);
dist_fit=dist(200:end);
modelstr='a*x^(0.5)+b*x^(1.5)+c';
model = fittype (modelstr);
Value = 'LinearLeastSquares';
opts = fitoptions('method', Value);
[fresult,gof,output] = fit(t_fit',dist_fit',model,opts);
if plotme>0
    figure('name','log scale')
    plot(t_fit,dist_fit,'ko-','LineWidth',2)
    hold on
    plot(fresult,'g')
    hold off
    set(gca,'FontSize',16)
    xlabel('t','FontSize',16)
    ylabel('\delta','FontSize',16)
    %--------------------------------------
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
    dist=dist-min(dist);
    dist2=dist(20:end);
    loglog(t2,dist2)
end
%---------Function return values ------------
a_fit=fresult.a
b_fit=fresult.b
c_fit=fresult.c
R_fit=gof.rsquare




