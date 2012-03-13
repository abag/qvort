function mesh_urms_all(start,finish,skip,smooth_fac,option)
if nargin<3
  disp('I do not have enough arguements to run, please supply start, final, skip')
  return
elseif nargin==3
  smooth_fac=1;
  option='3d';
elseif nargin==4
  option='3d';
elseif nargin>5
  disp('You have given me too many arguements, please consult the script')
  return 
end
counter=1;
for i=start:skip:finish
  switch option
  case '3d'
    [t urms nurms]=mesh_urms(i);
    tt(i)=t;
    u2(i)=urms;
    nu2(i)=nurms;
  case '2d'
    [urms2 nurms2]=two_dim_urms(i);
    u2(i)=urms2;
    nu2(i)=nurms2;
  end
  counter=counter+1;
end
disp('saving to u2.mat file')
save u2.mat u2
switch option
  case '3d'
    figure('Name','urms calculated from 3d mesh')
    plot(tt,smooth(u2,smooth_fac),'k','LineWidth',2)
  case '2d'
    figure('Name','urms calculated from 2d mesh')
    plot(smooth(u2,smooth_fac),'k','LineWidth',2)
end
set(gca,'FontSize',14)
xlabel('t','FontSize',14)
ylabel('E','FontSize',14)
