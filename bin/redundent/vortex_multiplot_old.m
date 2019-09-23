%secondary  vortex plotting routine, run this with a filenumber as an input, e.g.
%vortex_plot(1)
%the routine will also take the follwing options in the form vortex_plot(n,'option1','option2') e.g.
%options are:
%            rough: plots particles only advised for a large number of particles
%            annotate: add information on box size and time to plot
%            thin_line: if the tangle is very dense use this
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

%set options based on varargin
rough=0 ;  printit=0  ; annotate=0 ; thin_line=0 ; ploteps=0 ;
for i=1:optargin
   dummy_arg=cell2str(varargin(i));
   [dummy_arg value]=strtok(dummy_arg);
  switch  dummy_arg
    case 'rough'
      rough=1;
    case 'eps'
      ploteps=1;
    case 'print'
      printit=1;
      disp('printing to file')
      figure('visible', 'off')
    case 'annotate'
      annotate=1;
    case 'thin_line'
      thin_line=1;  
      otherwise
      disp('invalid option in input arguements')
      disp('aborting code and printing help:')
      help vortex_plot
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
  subplot(2,2,1)
    title('tangle')
    plot3(x,y,z,'.')
    axis([-dims(2)/2 dims(2)/2 -dims(2)/2 dims(2)/2 -dims(2)/2 dims(2)/2]);
    box on
  subplot(2,2,2)
    title('xy-plane')
    plot3(x,y,z,'.')
    axis([-dims(2)/2 dims(2)/2 -dims(2)/2 dims(2)/2 -dims(2)/2 dims(2)/2]);
    view(90,0)
    box on
  subplot(2,2,3)
    title('yz-plane')
    plot3(x,y,z,'.')
    axis([-dims(2)/2 dims(2)/2 -dims(2)/2 dims(2)/2 -dims(2)/2 dims(2)/2]);
    view(0,90)
    box on
  subplot(2,2,4)
    title('yz-plane')
    plot3(x,y,z,'.')
    axis([-dims(2)/2 dims(2)/2 -dims(2)/2 dims(2)/2 -dims(2)/2 dims(2)/2]);
    view(200,30)
    box on
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
    if (dist<4.*dims(1))
      if thin_line==1
        subplot(2,2,1)
        plot3(dummy_x(1:2,1),dummy_x(1:2,2),dummy_x(1:2,3),'-k','LineWidth',.5)
        hold on
        subplot(2,2,2)
        plot3(dummy_x(1:2,1),dummy_x(1:2,2),dummy_x(1:2,3),'-k','LineWidth',.5)
        hold on
        subplot(2,2,3)
        plot3(dummy_x(1:2,1),dummy_x(1:2,2),dummy_x(1:2,3),'-k','LineWidth',.5)                
        hold on
        subplot(2,2,4)
        plot3(dummy_x(1:2,1),dummy_x(1:2,2),dummy_x(1:2,3),'-k','LineWidth',.5)  
        hold on
      else
        subplot(2,2,1)
        plot3(dummy_x(1:2,1),dummy_x(1:2,2),dummy_x(1:2,3),'-k','LineWidth',2)
        hold on
        subplot(2,2,2)
        plot3(dummy_x(1:2,1),dummy_x(1:2,2),dummy_x(1:2,3),'-k','LineWidth',2)
        hold on
        subplot(2,2,3)
        plot3(dummy_x(1:2,1),dummy_x(1:2,2),dummy_x(1:2,3),'-k','LineWidth',2)
        hold on
        subplot(2,2,4)
        plot3(dummy_x(1:2,1),dummy_x(1:2,2),dummy_x(1:2,3),'-k','LineWidth',2)
        hold on
      end  
    end
  end
end
for i=1:4
  subplot(2,2,i)
  axis([-dims(2)/2 dims(2)/2 -dims(2)/(2*dims(7)) dims(2)/(2*dims(7)) -dims(2)/(2*dims(7)) dims(2)/(2*dims(7))]);
  daspect([1 dims(7) dims(7)])
  box on
  hold off
  if annotate==1
    if i==2
      s1='t=';
      s2=num2str(time);
      str=strcat(s1,s2);
      text(-.6*dims(2),.6*dims(2),.7*dims(2),str,'FontSize',14)
    end
  else
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    set(gca,'ztick',[])
  end
  rotate3d on
  set(gca,'FontSize',14)
end
subplot(2,2,1);title('tangle 1')
subplot(2,2,2);view(90,0);title('xz-plane')
subplot(2,2,3);view(0,90);title('xy-plane')
subplot(2,2,4);view(200,30);title('tangle 2')
if printit==1
  if ploteps==1 
    disp(sprintf('printing to vortex_out%04d.eps',filenumber))
    fOUT=sprintf('vortex_out%04d.eps',filenumber);
    print('-deps', fOUT)
  else
    disp(sprintf('printing to vortex_out%04d.png',filenumber))
    fOUT=sprintf('vortex_out%04d.png',filenumber);
    print('-dpng', fOUT) 
  end
end
