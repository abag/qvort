function mesh_urms_all(start,finish,skip,option)
if nargin<3
  disp('I do not have enough arguements to run, please supply start, final, skip')
  return
else if nargin==3
  option='3d'
else
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
    u22(i)=urms2;
    nu22(i)=nurms2;
  end
  counter=counter+1;
end
switch option
  case '3d'
    figure('Name','urms calculated from 3d mesh')
    plot(tt,u2,'k','LineWidth',2)
  case '2d'
    figure('Name','urms calculated from 2d mesh')
    plot(u22,'k','LineWidth',2)
end
set(gca,'FontSize',14)
xlabel('t','FontSize',14)
ylabel('E','FontSize',14)
