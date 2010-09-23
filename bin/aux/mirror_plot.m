function mirror_plot
filename='./data/mirror.dat';
%we set the dimensions of the box here
%this is overridden if we have periodic B.C.'s
dims=load('./data/dims.log');

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

%now create vectors to plot%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j=1:number_of_particles
  if round(f(j))==0
  else
    dummy_x(1,1)=x(j);
    dummy_x(2,1)=x(round(f(j)));
    dummy_x(1,2)=y(j);
    dummy_x(2,2)=y(round(f(j)));
    dummy_x(1,3)=z(j);
    dummy_x(2,3)=z(round(f(j)));
    pos=dummy_x(1,:);
    arrow=4*(dummy_x(2,:)-dummy_x(1,:));
    %plot3(dummy_x(1:2,1),dummy_x(1:2,2),dummy_x(1:2,3),'-m','LineWidth',2.0)
    arrow3d(pos,arrow,'k', 0.9)
    hold on
  end
end
boxx=[0 1 1 0 0 0;1 1 0 0 1 1;1 1 0 0 1 1;0 1 1 0 0 0]*dims(2)-dims(2)/2.;
boxy=[0 0 1 1 0 0;0 1 1 0 0 0;0 1 1 0 1 1;0 0 1 1 1 1]*dims(2)-dims(2)/2.;
boxz=[0 0 0 0 0 1;0 0 0 0 0 1;1 1 1 1 0 1;1 1 1 1 0 1]*dims(2)-dims(2)/2.;
for i=1:6
    h=patch(boxx(:,i),boxy(:,i),boxz(:,i),'b');
    alpha(h,0.1)
    set(h,'edgecolor','w');
end
axis equal