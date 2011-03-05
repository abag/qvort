function quasi_trace(varargin)
optargin = size(varargin,2);
disp('plotting a trace of quasi-particles paths')
load ./data/dims.log
pnumber=dims(6);
disp(sprintf('there are %04d quasi particles in the code',pnumber))
ts=load('data/ts.log');
filenumber=length(ts);
disp(sprintf('plotting trajectories from %04d snapshots',filenumber))
surface=0; dark=0; printit=0; 
for i=1:optargin
  switch cell2str(varargin(i))
    case 'surface'
      surface=1;
    case 'dark'
      dark=1;
    case 'print'
      printit=1;
      disp('printing to file')
      figure('visible', 'off')
    otherwise
      disp('invalid option in input arguements')
      disp('aborting code and printing help:')
      help particle_trace
      return
  end
end
%get the dimensions information from dims.log
dims=load('./data/dims.log');
c=colormap(jet(pnumber));
for i=1:filenumber
  filename=sprintf('data/quasi_par%04d.log',i);
  fid=fopen(filename);
  if dims(4)==1
    %binary read
    t=fread(fid,1,'float64');
    number_of_particles=fread(fid,1,'int');
    x=fread(fid,number_of_particles,'float64');
    y=fread(fid,number_of_particles,'float64');
    z=fread(fid,number_of_particles,'float64');
    for j=1:number_of_particles
      parray(j,i,1)=x(j);
      parray(j,i,2)=y(j);
      parray(j,i,3)=z(j);
    end
  else 
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
      parray(j,i,1)=x(j);
      parray(j,i,2)=y(j);
      parray(j,i,3)=z(j);
    end
  end
  fclose(fid);
end
for j=1:pnumber
    for i=1:filenumber-1
      dist=sqrt((parray(j,i,1)-parray(j,i+1,1))^2+(parray(j,i,2)-parray(j,i+1,2))^2+(parray(j,i,3)-parray(j,i+1,3))^2);
      if (dist<dims(2)/2.) || (dims(2)==0.)
        if surface==1
          [x1 y1 z1]=cylind(dims(1)/3,20, squeeze(parray(j,i,:))',squeeze(parray(j,i+1,:))');
          h=surf(x1,y1,z1);
          set(h,'FaceColor',c(j,:),'EdgeColor',c(j,:),'FaceAlpha',0.5,'EdgeAlpha',0.1) ;
        else
          plot3([parray(j,i,1) parray(j,i+1,1)],[parray(j,i,2) parray(j,i+1,2)],[parray(j,i,3) parray(j,i+1,3)],'Color',c(j,:),'LineWidth',2)
        end 
      end
        hold on
    end
end
%get the dimensions information from dims.log
if (dims(2)>0.)
  axis([-dims(2)/2 dims(2)/2 -dims(2)/2 dims(2)/2 -dims(2)/2 dims(2)/2]); 
  box on
else
  axis([-0.1 0.1 -0.1 0.1 -0.1 0.1]);      
end
if printit==1
  disp('printing to particle_trace.png')
  print('-dpng', 'particle_trace.png') 
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
set(gca,'FontSize',14)
rotate3d on
