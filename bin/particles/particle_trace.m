function vortex_plot(filenumber,pnumber)
disp(sprintf('plotting a trace of %04d particles trajectory',pnumber))
if nargin==1
  pnumber=1;
end
c=colormap(jet(pnumber));
for i=1:filenumber
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
    %if j<=pnumber
    %    plot3(x(j),y(j),z(j),'o','MarkerFaceColor',c(j,:))
    %    hold on
    %end
  end
end
dims=load('./data/dims.log');
for j=1:pnumber
    for i=1:filenumber-1
      dist=sqrt((parray(j,i,1)-parray(j,i+1,1))^2+(parray(j,i,2)-parray(j,i+1,2))^2+(parray(j,i,3)-parray(j,i+1,3))^2);
      if (dist<dims(2)/2.) || (dims(2)==0.)
        plot3([parray(j,i,1) parray(j,i+1,1)],[parray(j,i,2) parray(j,i+1,2)],[parray(j,i,3) parray(j,i+1,3)],'Color',c(j,:),'LineWidth',1.5)
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
