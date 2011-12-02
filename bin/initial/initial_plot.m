function initial_plot(varargin)
optargin = size(varargin,2);
%check filenumber has been set

filename='./data/var_initial.log';
%set options based on varargin
rough=0 ; linetrue=1 ;  dark=0 ; printit=0 ; overhead=0 ; show_points=0 ; overhead_xz=0;
annotate=0 ; thin_line=0 ; ploteps=0 ;
for i=1:optargin
  dummy_arg=cell2str(varargin(i));
  [dummy_arg value]=strtok(dummy_arg);
  switch  dummy_arg
    case 'rough'
      rough=1;
    case 'overhead'
      overhead=1;
    case 'overhead_xz'
      overhead_xz=1;
    case 'eps'
      ploteps=1;
    case 'print'
      printit=1;
      disp('printing to file')
      figure('visible', 'off')
    case 'show_points'
      show_points=1;
    case 'annotate'
      annotate=1
    case 'dark'
      dark=1;
    case 'smooth'
      linetrue=0;
    case 'thin_line'
      thin_line=1;
    otherwise
      disp('invalid option in input arguements')
      disp('aborting code')
      return
  end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%get the dimensions information from dims.log
dims=load('./data/dims.log');
fid=fopen(filename);
if fid<0
  disp('var file does not exist, exiting script')
  return
end
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
if (rough==1)
  plot3(x,y,z,'.')
  if (dims(2)>0.)
    if overhead==1
      axis([-dims(2)/2 dims(2)/2 -dims(2)/2 dims(2)/2]);
    else
      axis([-dims(2)/2 dims(2)/2 -dims(2)/2 dims(2)/2 -dims(2)/2 dims(2)/2]);
      box on
    end
  else
    if overhead==1
      axis([-box_size box_size -box_size box_size]);
    else
      axis([-box_size box_size -box_size box_size -box_size box_size]);
      box on
    end
  end
  rotate3d on
  return
end
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
    if (dist<0.5*dims(2))
      if linetrue==1
        if dark==1
          
          if show_points==1
            if (thin_line==1)
              plot3(dummy_x(1:2,1),dummy_x(1:2,2),dummy_x(1:2,3),'-om','LineWidth',.5,'MarkerFaceColor','m')
            else
              plot3(dummy_x(1:2,1),dummy_x(1:2,2),dummy_x(1:2,3),'-om','LineWidth',2.0)
            end
          else
            if (thin_line==1)
              plot3(dummy_x(1:2,1),dummy_x(1:2,2),dummy_x(1:2,3),'-m','LineWidth',.5)
            else
              plot3(dummy_x(1:2,1),dummy_x(1:2,2),dummy_x(1:2,3),'-m','LineWidth',2.0)
            end
          end
          
        else
          
          if show_points==1
            if thin_line==1
              plot3(dummy_x(1:2,1),dummy_x(1:2,2),dummy_x(1:2,3),'-ok','LineWidth',.5,'MarkerFaceColor','k')
            else
              plot3(dummy_x(1:2,1),dummy_x(1:2,2),dummy_x(1:2,3),'-ok','LineWidth',2)
            end
          else
            if thin_line==1
              plot3(dummy_x(1:2,1),dummy_x(1:2,2),dummy_x(1:2,3),'-k','LineWidth',.1)
            else
              plot3(dummy_x(1:2,1),dummy_x(1:2,2),dummy_x(1:2,3),'-k','LineWidth',2)
            end
          end
        end
      else
        [x1 y1 z1]=cylind(dims(1)/2,20, dummy_x(1,1:3),dummy_x(2,1:3));
        h=surf(x1,y1,z1);
        if dark==1
          
          set(h,'FaceColor','m','EdgeColor','m','FaceAlpha',0.5,'EdgeAlpha',0.1) ;
          
        else
          
          if (j==-1) %set this off by default
            %pick out a particle in a particular colour?
            set(h,'FaceColor','k','EdgeColor','b','FaceAlpha',0.5,'EdgeAlpha',0.1) ;
          else
            set(h,'FaceColor','k','EdgeColor','k','FaceAlpha',0.2,'EdgeAlpha',0.1) ;
          end
          
        end
      end
      hold on
    end
  end
end
if (dims(2)>0.)
  if overhead==1
    axis([-dims(2)/2 dims(2)/2 -dims(2)/(2*dims(7)) dims(2)/(2*dims(7))]);
  elseif overhead_xz==1
    axis([-dims(2)/2 dims(2)/2 -dims(2)/(2*dims(7)) dims(2)/(2*dims(7)) -dims(2)/(2*dims(7)) dims(2)/(2*dims(7))]);
    view(0,0)
  else
    axis([-dims(2)/2 dims(2)/2 -dims(2)/(2*dims(7)) dims(2)/(2*dims(7)) -dims(2)/(2*dims(7)) dims(2)/(2*dims(7))]);
    daspect([1 dims(7) dims(7)])
    box on
  end
else
  if overhead==1
    axis([-box_size box_size -box_size box_size]);
  elseif overhead_xz==1
    axis([-dims(2)/2 dims(2)/2 -dims(2)/(2*dims(7)) dims(2)/(2*dims(7)) -dims(2)/(2*dims(7)) dims(2)/(2*dims(7))]);
    view(0,0)
  else
    axis([-dims(2)/2 dims(2)/2 -dims(2)/(2*dims(7)) dims(2)/(2*dims(7)) -dims(2)/(2*dims(7)) dims(2)/(2*dims(7))]);
    daspect([1 dims(7) dims(7)])
    box on
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
if annotate==1
  s1='t=';
  s2=num2str(time);
  str=strcat(s1,s2);
  if (dims(2)>0.)
    text(-.6*dims(2),.6*dims(2),.7*dims(2),str,'FontSize',14)
  else
    text(-0.06,0.06,0.07,str,'FontSize',16)
  end
end
if printit==1
  if ploteps==1
    disp(sprintf('printing to vortex_out%04d.eps',filenumber))
    fOUT=sprintf('vortex_out%04d.eps',filenumber)
    print('-deps', fOUT)
  else
    disp(sprintf('printing to vortex_out%04d.png',filenumber))
    fOUT=sprintf('vortex_out%04d.png',filenumber)
    print('-dpng', fOUT)
  end
end
if annotate==0
  set(gca,'xtick',[])
  set(gca,'ytick',[])
  set(gca,'ztick',[])
end
rotate3d on
set(gca,'FontSize',14)
%text(1.,10.5,0.55,str)
%view(-18,17)

