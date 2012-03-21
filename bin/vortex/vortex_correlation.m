function vortex_correlation(filenumber)
%check filenumber has been set
if exist('filenumber')==0
  disp('you have not set filnumber')
  disp('aborting code and type "help vortex_plot" for more options')
  return
end
filename=sprintf('data/var%04d.log',filenumber);
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
%loaded in data now loop over all points
counter=1;
for j=1:100 %number_of_particles
    j/100
    for i=j:100 %number_of_particles
        if round(f(j))==0
            continue
        end
         if round(f(i))==0
            continue
        end
        dummy_x(1,1)=x(j);
        dummy_x(2,1)=x(round(f(j)));
        dummy_x(1,2)=y(j);
        dummy_x(2,2)=y(round(f(j)));
        dummy_x(1,3)=z(j);
        dummy_x(2,3)=z(round(f(j)));
        dummy_y(1,1)=x(i);
        dummy_y(2,1)=x(round(f(i)));
        dummy_y(1,2)=y(i);
        dummy_y(2,2)=y(round(f(i)));
        dummy_y(1,3)=z(i);
        dummy_y(2,3)=z(round(f(i)));
        distx=sqrt((dummy_x(1,1)-dummy_x(2,1))^2+(dummy_x(1,2)-dummy_x(2,2))^2+(dummy_x(1,3)-dummy_x(2,3))^2);
        if (distx>0.5*dims(2))
          continue
        end
        disty=sqrt((dummy_y(1,1)-dummy_y(2,1))^2+(dummy_y(1,2)-dummy_y(2,2))^2+(dummy_y(1,3)-dummy_y(2,3))^2);
        if (disty>0.5*dims(2))
          continue
        end
        dist=norm(0.5*(dummy_x(1,:)+dummy_x(2,:))-0.5*(dummy_y(1,:)+dummy_y(2,:)));
        uni_tangx=(dummy_x(2,:)-dummy_x(1,:))/distx;
        uni_tangy=(dummy_y(2,:)-dummy_y(1,:))/disty;
        Qxx(counter,1)=dist;
        Qxx(counter,2)=uni_tangx(1)*uni_tangy(1);
        counter=counter+1;
    end
end
min_r=min(Qxx(:,1));
max_r=max(Qxx(:,1));
nbins=10;
Qxx_avg(1:nbins-1)=0.;
counter2(1:nbins-1)=0.;
for i=1:length(Qxx)
    for j=1:nbins-1
        if (Qxx(i,1)>(((j-1)*(max_r-min_r)/(nbins-1))+min_r)) && (Qxx(i,1)<(((j)*(max_r-min_r)/(nbins-1))+min_r))
            Qxx_avg(j)=Qxx_avg(j)+Qxx(i,2);
            counter2(j)=counter2(j)+1;
        end
    end
end
for i=1:nbins-1
    if (counter2(i)>0.)
        Qxx_avg(i)=Qxx_avg(i)/counter2(i);
    else
        Qxx_avg(i)=0.;
    end
end
r=(1:(nbins-1))*sqrt(3*dims(2)^2)/nbins;
plot(r,Qxx_avg)
