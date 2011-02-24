function particle_shadow(filenumber)
disp('plotting 2 particle trajecties with shadow effect')
if nargin==1
  pnumber=1;
end
if filenumber<2
  disp('setting filenumber to be 2 to achieve shadow effect')
  filenumber=2
end
start_fnumber=max(1,filenumber-20);
for i=start_fnumber:filenumber
  filename=sprintf('data/par%04d.log',i);
  fid=fopen(filename);
  %read the time
  tline=fgetl(fid);
  dummy=textscan(tline, '%f');
  time=dummy{:};
  %how many particles
  tline=fgetl(fid);
  dummy=textscan(tline, '%d');
  number_of_particles=dummy{:};
  if number_of_particles<2
    disp('there is only one particle - exiting script')
    return
  end
  %get the particles%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  %particle 1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  tline=fgetl(fid);
  dummy=textscan(tline, '%f');
  dummy_vect=dummy{:};
  x=dummy_vect(1);
  y=dummy_vect(2);
  z=dummy_vect(3);
  p1(i,1)=x;
  p1(i,2)=y;
  p1(i,3)=z;
  %particle 1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  tline=fgetl(fid);
  dummy=textscan(tline, '%f');
  dummy_vect=dummy{:};
  x=dummy_vect(1);
  y=dummy_vect(2);
  z=dummy_vect(3);
  p2(i,1)=x;
  p2(i,2)=y;
  p2(i,3)=z;
end
%get the dimensions information from dims.log
dims=load('./data/dims.log');
%plot the first particle
s=size(p1);
for i=1:s(1)-1
  alpha_lev=(i/(s(1)-1))^10;
  dist=sqrt((p1(i,1)-p1(i+1,1))^2+(p1(i,2)-p1(i+1,2))^2+(p1(i,3)-p1(i+1,3))^2);
  if (dist<dims(2)/2.) || (dims(2)==0.)
    [x1 y1 z1]=cylind(0.00045,20, p1(i,1:3),p1(i+1,1:3));
    h=surf(x1,y1,z1);
    set(h,'FaceColor','k','EdgeColor','k','FaceAlpha',alpha_lev,'EdgeAlpha',alpha_lev) ;
  end
  %plot3(p1(i:i+1,1),p1(i:i+1,2),p1(i:i+1,3),'LineWidth',width)
  hold on
end
lighting phong
camlight
if (dims(2)>0.)
  axis([-dims(2)/2 dims(2)/2 -dims(2)/2 dims(2)/2 -dims(2)/2 dims(2)/2]); 
  box on
else
  axis([-0.1 0.1 -0.1 0.1 -0.1 0.1]);      
end
