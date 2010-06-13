function vortex_plot(filenumber)
filename=sprintf('data/var%03d.log',filenumber);
%some options
linetrue=0; %if 1 plots a line, else plots a thin cylinder
dark=1; %if 1 plots in dark
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fid=fopen(filename);
%read the time
tline=fgetl(fid);
dummy=textscan(tline, '%f');
time=dummy{:};
%how many particles
tline=fgetl(fid);
dummy=textscan(tline, '%d');
number_of_particles=dummy{:};

%get the particles%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j=1:number_of_particles
  tline=fgetl(fid);
  dummy=textscan(tline, '%f');
  dummy_vect=dummy{:};
  x(j)=dummy_vect(1);
  y(j)=dummy_vect(2);
  z(j)=dummy_vect(3);
  f(j)=dummy_vect(4);
end
f=uint16(f);
%get the dimensions information from dims.log
dims=load('./data/dims.log');

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
    dist=sqrt((dummy_x(1,1)-dummy_x(2,1))^2+(dummy_x(1,2)-dummy_x(2,2))^2+(dummy_x(1,3)-dummy_x(2,3))^2);
    if (dist<3.*dims(1))
      if linetrue==1
        if dark==1
          plot3(dummy_x(1:2,1),dummy_x(1:2,2),dummy_x(1:2,3),'-m','LineWidth',2.0)
        else
          plot3(dummy_x(1:2,1),dummy_x(1:2,2),dummy_x(1:2,3),'-k','LineWidth',2.0)
        end
      else
        [x1 y1 z1]=cylind(0.00045,20, dummy_x(1,1:3),dummy_x(2,1:3));
        h=surf(x1,y1,z1);
        if dark==1
          set(h,'FaceColor','m','EdgeColor','m','FaceAlpha',0.5,'EdgeAlpha',0.1) ;
        else
          set(h,'FaceColor','k','EdgeColor','k','FaceAlpha',0.5,'EdgeAlpha',0.1) ;
        end
      end
      if (dims(2)>0.)
        axis([-dims(2)/2 dims(2)/2 -dims(2)/2 dims(2)/2 -dims(2)/2 dims(2)/2]); 
        box on
      else
        axis([-0.1 0.1 -0.1 0.1 -0.1 0.1]);
      end
      hold on
    end
  end
end
if dark==1
  whitebg('k')
  set(gcf,'InvertHardcopy','off');
else
  whitebg('w')
end
lighting phong
camlight
hold off
grid on
s1='t=';
s2=num2str(time);
str=strcat(s1,s2);
%text(-0.06,0.06,0.07,str,'FontSize',16)
%text(1.,10.5,0.55,str)
%view(-18,17)

