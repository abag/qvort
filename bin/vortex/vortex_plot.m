%the main vortex plotting routine, run this with a filenumber as an input, e.g.
%vortex_plot(1)
%the routine will also take the follwing options in the form vortex_plot(n,'option1','option2') e.g.
%options are:
%            rough: plots particles only advised for a large number of particles
%            line: looks better than above a nice compromise
%            dark: add a night time theme!
%            annotate: add information on box size and time to plot
%	     rainbow: colour code the vortex according to velocity/magnetic field
%	     rainbow: colour code the vortex according to log velocity/magnetic field
%	     magnetic: colour code a magnetic flux tube according to B, will also affext the width of tubes
%            vmax= value:maximum field strength
%            vmin= value:maximum field strength
%            overhead: angle the plot overhead
%            overhead_xz: angle the plot overhead xz plane
%            show_points: will only work if line is set, shows points as well as lines, ignored if rainbow set
%            print: print to file rather than screen
%            eps: if print is set and eps is set then output eps files
%
function vortex_plot(filenumber,varargin)
optargin = size(varargin,2);
%check filenumber has been set
if exist('filenumber')==0
  disp('you have not set filnumber')
  disp('aborting code and type "help vortex_plot" for more options')
  return
end
filename=sprintf('data/var%04d.log',filenumber);
%we set the dimensions of the box here
%this is overridden if we have periodic B.C.'s
box_size=.005 ;
%set options based on varargin
rough=0 ; linetrue=0 ; rainbow=0 ; dark=0 ; printit=0 ; overhead=0  ; magnetic=0; show_points=0 ; sph_associated=0 ; overhead_xz=0;
log_rainbow=0 ; annotate=0 ; thin_line=0 ; ploteps=0 ;
%empty the vmax/min values
v_max=[] ; v_min=[] ;
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
    case 'rainbow'
      rainbow=1;
    case 'log_rainbow'
      log_rainbow=1;
    case 'show_points'
      show_points=1;
    case 'vmax='
       v_max=str2num(value);
     case 'vmin='
       v_min=str2num(value);
    case 'sph_associated'
      sph_associated=1;
    case 'magnetic'
      rainbow=0; %switchoff-rainbow
      magnetic=1;
    case 'annotate'
      annotate=1
    case 'dark'
      dark=1;
    case 'line'
      linetrue=1;
    case 'thin_line'
      linetrue=1;
      thin_line=1;  
      otherwise
      disp('invalid option in input arguements')
      disp('aborting code and printing help:')
      help vortex_plot
      return
  end
end
if (magnetic==1) 
   if isempty(v_max)==1
       v_max=1.;
   end
    if isempty(v_min)==1
       v_min=0.;
    end
    if (v_min>v_max)
      disp('you have set v_min>v_max ; exiting script')
      return
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%get the dimensions information from dims.log
dims=load('./data/dims.log');
if dims(4)==1
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
  u=fread(fid,number_of_particles,'float64');
  u2=fread(fid,number_of_particles,'float64');
else 
  fid=fopen(filename);
  if fid<0
      disp('var file does not exist, exiting script')
      return
  end
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
    u(j)=dummy_vect(5);
    u2(j)=dummy_vect(6);
  end
  f=uint16(f);
end
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
if rainbow==1
  %scale velocity into a colormap
  store_caxis=([min(u(u>0)) max(u)]);
  u=u-min(u(u>0));
  rainbow_scale=199/max(u) ;
  u=u*rainbow_scale;
  rainbowcmap=colormap(jet(200));
elseif log_rainbow==1
  rainbow=1 ;
  %scale velocity into a colormap
  store_caxis=([min(log10(u(u>0))) max(log10(u))]);
  u=log10(u)-min(log10(u(u>0)));
  rainbow_scale=199/max(u) ;
  u=u*rainbow_scale;
  rainbowcmap=colormap(jet(200));
elseif magnetic==1
  if (max(u)>v_max)
      disp('v_max is too small ; exiting script')
      return
  end 
  if (min(u)<v_min)
      disp('v_min is too large ; exiting script')
      return
  end 
  %scale field into a colormap
  twid=1./u;
  store_caxis=([v_min v_max]);
  u=u-v_min;
  rainbow_scale=299/v_max ;
  u=u*rainbow_scale;
  rainbowcmap=colormap(jet(300));
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
    if (dist<4.*dims(1))
      if linetrue==1
        if dark==1
          if rainbow==1
            if u(j)==0
              u(j)=1;
            end
            if (thin_line==1)
              plot3(dummy_x(1:2,1),dummy_x(1:2,2),dummy_x(1:2,3),'-','Color',rainbowcmap(max(1,ceil(u(j))),:),'LineWidth',.5)
            else
              plot3(dummy_x(1:2,1),dummy_x(1:2,2),dummy_x(1:2,3),'-','Color',rainbowcmap(max(1,ceil(u(j))),:),'LineWidth',2.0)
            end
          elseif magnetic==1
            if u(j)==0
              u(j)=1;
            end
            if (thin_line==1)
              plot3(dummy_x(1:2,1),dummy_x(1:2,2),dummy_x(1:2,3),'-','Color',rainbowcmap(max(1,ceil(u(j))),:),'LineWidth',.5)
            else
              plot3(dummy_x(1:2,1),dummy_x(1:2,2),dummy_x(1:2,3),'-','Color',rainbowcmap(max(1,ceil(u(j))),:),'LineWidth',2.0)
            end
          else
            if show_points==1
              if (thin_line==1)
                plot3(dummy_x(1:2,1),dummy_x(1:2,2),dummy_x(1:2,3),'-om','LineWidth',.5)
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
          end
        else
          if rainbow==1
            if u(j)==0
              u(j)=1;
            end
            if (thin_line==1)  
              plot3(dummy_x(1:2,1),dummy_x(1:2,2),dummy_x(1:2,3),'-','Color',rainbowcmap(max(1,ceil(u(j))),:),'LineWidth',.5)
            else
              plot3(dummy_x(1:2,1),dummy_x(1:2,2),dummy_x(1:2,3),'-','Color',rainbowcmap(max(1,ceil(u(j))),:),'LineWidth',2.0)   
            end
          elseif magnetic==1
            if u(j)==0
              u(j)=1;
            end
            if (thin_line==1)  
              plot3(dummy_x(1:2,1),dummy_x(1:2,2),dummy_x(1:2,3),'-','Color',rainbowcmap(max(1,ceil(u(j))),:),'LineWidth',.5)
            else
              plot3(dummy_x(1:2,1),dummy_x(1:2,2),dummy_x(1:2,3),'-','Color',rainbowcmap(max(1,ceil(u(j))),:),'LineWidth',2.0)   
            end
          else
            if show_points==1
              if thin_line==1
                plot3(dummy_x(1:2,1),dummy_x(1:2,2),dummy_x(1:2,3),'-ok','LineWidth',.5)
              else
                plot3(dummy_x(1:2,1),dummy_x(1:2,2),dummy_x(1:2,3),'-ok','LineWidth',2)   
              end
            else
              if sph_associated==1
                if u2(j)<1
                  plot3(dummy_x(1:2,1),dummy_x(1:2,2),dummy_x(1:2,3),'-k','LineWidth',2)
                else
                  plot3(dummy_x(1:2,1),dummy_x(1:2,2),dummy_x(1:2,3),'-m','LineWidth',2)  
                end
              else
                if thin_line==1
                  plot3(dummy_x(1:2,1),dummy_x(1:2,2),dummy_x(1:2,3),'-k','LineWidth',.1)  
                else
                  plot3(dummy_x(1:2,1),dummy_x(1:2,2),dummy_x(1:2,3),'-k','LineWidth',2)  
                end
              end
            end 
          end
        end
      else
        if magnetic==1
          [x1 y1 z1]=cylind(twid(j)*dims(1)/5,20, dummy_x(1,1:3),dummy_x(2,1:   3));
        else
            [x1 y1 z1]=cylind(dims(1)/8,20, dummy_x(1,1:3),dummy_x(2,1:3));
        end
        h=surf(x1,y1,z1);
        if dark==1
          if rainbow==1
            if u(j)==0
                u(j)=1;
            end
            set(h,'FaceColor',rainbowcmap(max(1,ceil(u(j))),:),'EdgeColor',rainbowcmap(max(1,ceil(u(j))),:),'FaceAlpha',0.5,'EdgeAlpha',0.1) ;
          else
            set(h,'FaceColor','m','EdgeColor','m','FaceAlpha',0.5,'EdgeAlpha',0.1) ;
          end
        else
          if rainbow==1
            if u(j)==0
                u(j)=1;
            end
            set(h,'FaceColor',rainbowcmap(max(1,ceil(u(j))),:),'EdgeColor',rainbowcmap(max(1,ceil(u(j))),:),'FaceAlpha',0.5,'EdgeAlpha',0.1) ;
          elseif magnetic==1
            if u(j)==0
                u(j)=1;
            end
            set(h,'FaceColor',rainbowcmap(max(1,ceil(u(j))),:),'EdgeColor',rainbowcmap(max(1,ceil(u(j))),:),'FaceAlpha',0.5,'EdgeAlpha',0.1) ;
          else
            if (j==-1) %set this off by default
              %pick out a particle in a particular colour?
              set(h,'FaceColor','k','EdgeColor','b','FaceAlpha',0.5,'EdgeAlpha',0.1) ;
            else
              set(h,'FaceColor','k','EdgeColor','k','FaceAlpha',0.2,'EdgeAlpha',0.1) ; 
            end
          end
        end
      end
      hold on
    end
  end
end
if (dims(2)>0.)
  if overhead==1
    axis([-dims(2)/2 dims(2)/2 -dims(2)/(2*dims(9)) dims(2)/(2*dims(9))]);
  elseif overhead_xz==1
    axis([-dims(2)/2 dims(2)/2 -dims(2)/(2*dims(9)) dims(2)/(2*dims(9)) -dims(2)/(2*dims(9)) dims(2)/(2*dims(9))]);
    view(0,0)
  else
    axis([-dims(2)/2 dims(2)/2 -dims(2)/(2*dims(9)) dims(2)/(2*dims(9)) -dims(2)/(2*dims(9)) dims(2)/(2*dims(9))]);
    daspect([1 dims(9) dims(9)])
    box on
  end
else
  if overhead==1
    axis([-box_size box_size -box_size box_size]);
  elseif overhead_xz==1
    axis([-dims(2)/2 dims(2)/2 -dims(2)/(2*dims(9)) dims(2)/(2*dims(9)) -dims(2)/(2*dims(9)) dims(2)/(2*dims(9))]);
    view(0,0)
  else
    axis([-dims(2)/2 dims(2)/2 -dims(2)/(2*dims(9)) dims(2)/(2*dims(9)) -dims(2)/(2*dims(9)) dims(2)/(2*dims(9))]);
    daspect([1 dims(9) dims(9)])
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
if rainbow==1
  caxis(store_caxis)
  colorbar
end
if magnetic==1
  caxis(store_caxis)
  colorbar
end
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

