%how many distinct vortex filaments do we have in the simulation at
%present print these loops to separate files 
function vortex_kelvin_wave_analysis(filenumber)
%global dummy_pos
filename=sprintf('data/var%04d.log',filenumber);
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
zerocount=sum(f<=0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%find a starting point
for j=1:number_of_particles
  if f(j)~=0
    next=j;
    break
  end
end
%%%%%%%%%%%now do main loop%%%%%%%%%%%%%
line_count=0;
next_old=next;
counter(1:500)=0;
cmap=colormap(lines(500));
for l=1:500
  for i=1:number_of_particles
    line(l,i)=next;
    old_next=next;
    next=f(next);
    counter(l)=counter(l)+1;
    dummy_pos(counter(l),1)=x(old_next);
    dummy_pos(counter(l),2)=y(old_next);
    dummy_pos(counter(l),3)=z(old_next);
    if next==next_old
      break
      counter(l)
    end
    if next==0
      break
    end
  end
  dummy_pos(length(dummy_pos)+1,1:3)=dummy_pos(1,1:3);
  if counter(l)>15000
      disp('The length of the line found is over 15,000 particles, it is:')
      length(dummy_pos)
      temp_dummy_pos=dummy_pos;
      break
  end
  if sum(counter)<(number_of_particles-zerocount)
      for i=1:number_of_particles
          unique=true;
          for m=1:l
              for j=1:counter(m)
                  if i==line(m,j)
                      unique=false;
                  end
                  if f(i)==0
                      unique=false;
                  end
              end
          end
          if unique
              next=i;
              next_old=next;
          end
      end
  else
      break
  end
  clear dummy_pos
end
dummy_pos=temp_dummy_pos;
clear temp_dummy_pos
length(dummy_pos)
for k=1:(length(dummy_pos)-1)
    dist=sqrt((dummy_pos(k,1)-dummy_pos(k+1,1))^2+(dummy_pos(k,2)-dummy_pos(k+1,2))^2+(dummy_pos(k,3)-dummy_pos(k+1,3))^2);
    if dist>10*dims(1)
        disp('unfolding the loop')
        if dummy_pos(k,1)-dummy_pos(k+1,1)>10*dims(1)
            dummy_pos(k+1:length(dummy_pos),1)=dummy_pos(k+1:length(dummy_pos),1)+dims(2);
            disp('unfolding the loop in the x direction')
        end
        if dummy_pos(k,1)-dummy_pos(k+1,1)<-10*dims(1)
            dummy_pos(k+1:length(dummy_pos),1)=dummy_pos(k+1:length(dummy_pos),1)-dims(2);
            disp('unfolding the loop in the x direction')
        end
        if dummy_pos(k,2)-dummy_pos(k+1,2)>10*dims(1)
            dummy_pos(k+1:length(dummy_pos),2)=dummy_pos(k+1:length(dummy_pos),2)+dims(2);
            disp('unfolding the loop in the y direction')
        end
        if dummy_pos(k,2)-dummy_pos(k+1,2)<-10*dims(1)
            dummy_pos(k+1:length(dummy_pos),2)=dummy_pos(k+1:length(dummy_pos),2)-dims(2);
            disp('unfolding the loop in the y direction')
        end
        if dummy_pos(k,3)-dummy_pos(k+1,3)>10*dims(1)
            dummy_pos(k+1:length(dummy_pos),3)=dummy_pos(k+1:length(dummy_pos),3)+dims(2);
            disp('unfolding the loop in the +ve z direction')
        end
        if dummy_pos(k,3)-dummy_pos(k+1,3)<-10*dims(1)
            dummy_pos(k+1:length(dummy_pos),3)=dummy_pos(k+1:length(dummy_pos),3)-dims(2);
            disp('unfolding the loop in the -ve z direction')
        end
        plot3(dummy_pos(1:k,1),dummy_pos(1:k,2),dummy_pos(1:k,3),'Color','k','LineWidth',1.5);
    end
end
plot3(dummy_pos(:,1),dummy_pos(:,2),dummy_pos(:,3),'Color','k','LineWidth',1.5);
%%%%%%%%%%%%%%%%%%%CROSS VALIDATION%%%%%%%%%%%%%%%%%%%%%
%   num = 100;
%   spans = linspace(0.01,0.99,num);
%   sse = zeros(size(spans));
%   cp = cvpartition(100,'k',10);
%   for j=1:length(spans),
%       f = @(train,test) norm(test(:,2) -mylowess(train,test(:,1),spans(j)))^2;
%       sse(j) = sum(crossval(f,[1:length(dummy_pos(:,1)),dummy_pos(:,1)],'partition',cp));
%   end
%   [minsse,minj] = min(sse);
%   span = spans(minj);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
smooth_dummy_pos20(:,1)=smooth(1:length(dummy_pos(:,1)),dummy_pos(:,1),20,'moving');
smooth_dummy_pos20(:,2)=smooth(1:length(dummy_pos(:,2)),dummy_pos(:,2),20,'moving');
smooth_dummy_pos20(:,3)=smooth(1:length(dummy_pos(:,3)),dummy_pos(:,3),20,'moving');

smooth_dummy_pos5(:,1)=smooth(1:length(dummy_pos(:,1)),dummy_pos(:,1),5,'moving');
smooth_dummy_pos5(:,2)=smooth(1:length(dummy_pos(:,2)),dummy_pos(:,2),5,'moving');
smooth_dummy_pos5(:,3)=smooth(1:length(dummy_pos(:,3)),dummy_pos(:,3),5,'moving');

smooth_dummy_pos10(:,1)=smooth(1:length(dummy_pos(:,1)),dummy_pos(:,1),10,'moving');
smooth_dummy_pos10(:,2)=smooth(1:length(dummy_pos(:,2)),dummy_pos(:,2),10,'moving');
smooth_dummy_pos10(:,3)=smooth(1:length(dummy_pos(:,3)),dummy_pos(:,3),10,'moving');

smooth_dummy_pos15(:,1)=smooth(1:length(dummy_pos(:,1)),dummy_pos(:,1),15,'moving');
smooth_dummy_pos15(:,2)=smooth(1:length(dummy_pos(:,2)),dummy_pos(:,2),15,'moving');
smooth_dummy_pos15(:,3)=smooth(1:length(dummy_pos(:,3)),dummy_pos(:,3),15,'moving');

smooth_dummy_pos25(:,1)=smooth(1:length(dummy_pos(:,1)),dummy_pos(:,1),25,'moving');
smooth_dummy_pos25(:,2)=smooth(1:length(dummy_pos(:,2)),dummy_pos(:,2),25,'moving');
smooth_dummy_pos25(:,3)=smooth(1:length(dummy_pos(:,3)),dummy_pos(:,3),25,'moving');

smooth_dummy_pos30(:,1)=smooth(1:length(dummy_pos(:,1)),dummy_pos(:,1),30,'moving');
smooth_dummy_pos30(:,2)=smooth(1:length(dummy_pos(:,2)),dummy_pos(:,2),30,'moving');
smooth_dummy_pos30(:,3)=smooth(1:length(dummy_pos(:,3)),dummy_pos(:,3),30,'moving');

smooth_dummy_pos50(:,1)=smooth(1:length(dummy_pos(:,1)),dummy_pos(:,1),50,'moving');
smooth_dummy_pos50(:,2)=smooth(1:length(dummy_pos(:,2)),dummy_pos(:,2),50,'moving');
smooth_dummy_pos50(:,3)=smooth(1:length(dummy_pos(:,3)),dummy_pos(:,3),50,'moving');

smooth_dummy_pos100(:,1)=smooth(1:length(dummy_pos(:,1)),dummy_pos(:,1),100,'moving');
smooth_dummy_pos100(:,2)=smooth(1:length(dummy_pos(:,2)),dummy_pos(:,2),100,'moving');
smooth_dummy_pos100(:,3)=smooth(1:length(dummy_pos(:,3)),dummy_pos(:,3),100,'moving');

smooth_dummy_pos200(:,1)=smooth(1:length(dummy_pos(:,1)),dummy_pos(:,1),200,'moving');
smooth_dummy_pos200(:,2)=smooth(1:length(dummy_pos(:,2)),dummy_pos(:,2),200,'moving');
smooth_dummy_pos200(:,3)=smooth(1:length(dummy_pos(:,3)),dummy_pos(:,3),200,'moving');

hold on
plot3(smooth_dummy_pos20(:,1),smooth_dummy_pos20(:,2),smooth_dummy_pos20(:,3),'Color','r','LineWidth',1.5)
hold off
disp('the length of dummy_pos is:')
length(dummy_pos)
disp('the length of smooth_dummy_pos20 is:')
length(smooth_dummy_pos20)

for j=20:(length(dummy_pos)-20)
   A_20(j-19)=sqrt((dummy_pos(j,1)-smooth_dummy_pos20(j,1))^2+(dummy_pos(j,2)-smooth_dummy_pos20(j,2))^2+(dummy_pos(j,3)-smooth_dummy_pos20(j,3))^2);
   xi1_20(j-19)=sqrt((dummy_pos(j+1,1)-dummy_pos(j,1))^2+(dummy_pos(j+1,2)-dummy_pos(j,2))^2+(dummy_pos(j+1,3)-dummy_pos(j,3))^2);
   xi2_20(j-19)=sqrt((smooth_dummy_pos20(j+1,1)-smooth_dummy_pos20(j,1))^2+(smooth_dummy_pos20(j+1,2)-smooth_dummy_pos20(j,2))^2+(smooth_dummy_pos20(j+1,3)-smooth_dummy_pos20(j,3))^2);
end
for k=1:length(xi1_20)
    xi11_20(k)=sum(xi1_20(1:k));
end
xi1_interp20 = 0:dims(1)/2:xi11_20(end);
A_interp20 = interp1(xi11_20,A_20,xi1_interp20);
A_dash20 = gradient(A_interp20(10:end),dims(1)/2);
figure
plot(xi1_interp20(10:end),A_interp20(10:end))
disp('when the span is 20, the mean is:')
mean(abs(A_dash20))

for j=5:(length(dummy_pos)-5)
   A_5(j-4)=sqrt((dummy_pos(j,1)-smooth_dummy_pos5(j,1))^2+(dummy_pos(j,2)-smooth_dummy_pos5(j,2))^2+(dummy_pos(j,3)-smooth_dummy_pos5(j,3))^2);
   xi1_5(j-4)=sqrt((dummy_pos(j+1,1)-dummy_pos(j,1))^2+(dummy_pos(j+1,2)-dummy_pos(j,2))^2+(dummy_pos(j+1,3)-dummy_pos(j,3))^2);
   xi2_5(j-4)=sqrt((smooth_dummy_pos5(j+1,1)-smooth_dummy_pos5(j,1))^2+(smooth_dummy_pos5(j+1,2)-smooth_dummy_pos5(j,2))^2+(smooth_dummy_pos5(j+1,3)-smooth_dummy_pos5(j,3))^2);
end
for k=1:length(xi1_5)
    xi11_5(k)=sum(xi1_5(1:k));
end
xi1_interp5 = 0:dims(1)/2:xi11_5(end);
A_interp5 = interp1(xi11_5,A_5,xi1_interp5);
A_dash5 = gradient(A_interp5(10:end),dims(1)/2);
figure
plot(xi1_interp5(10:end),A_interp5(10:end))
disp('when the span is 5, the mean is:')
mean(abs(A_dash5))

for j=10:(length(dummy_pos)-10)
   A_10(j-9)=sqrt((dummy_pos(j,1)-smooth_dummy_pos10(j,1))^2+(dummy_pos(j,2)-smooth_dummy_pos10(j,2))^2+(dummy_pos(j,3)-smooth_dummy_pos10(j,3))^2);
   xi1_10(j-9)=sqrt((dummy_pos(j+1,1)-dummy_pos(j,1))^2+(dummy_pos(j+1,2)-dummy_pos(j,2))^2+(dummy_pos(j+1,3)-dummy_pos(j,3))^2);
   xi2_10(j-9)=sqrt((smooth_dummy_pos10(j+1,1)-smooth_dummy_pos10(j,1))^2+(smooth_dummy_pos10(j+1,2)-smooth_dummy_pos10(j,2))^2+(smooth_dummy_pos10(j+1,3)-smooth_dummy_pos10(j,3))^2);
end
for k=1:length(xi1_10)
    xi11_10(k)=sum(xi1_10(1:k));
end
xi1_interp10 = 0:dims(1)/2:xi11_10(end);
A_interp10 = interp1(xi11_10,A_10,xi1_interp10);
A_dash10 = gradient(A_interp10(10:end),dims(1)/2);
figure
plot(xi1_interp10(10:end),A_interp10(10:end))
disp('when the span is 10, the mean is:')
mean(abs(A_dash10))

for j=15:(length(dummy_pos)-15)
   A_15(j-14)=sqrt((dummy_pos(j,1)-smooth_dummy_pos15(j,1))^2+(dummy_pos(j,2)-smooth_dummy_pos15(j,2))^2+(dummy_pos(j,3)-smooth_dummy_pos15(j,3))^2);
   xi1_15(j-14)=sqrt((dummy_pos(j+1,1)-dummy_pos(j,1))^2+(dummy_pos(j+1,2)-dummy_pos(j,2))^2+(dummy_pos(j+1,3)-dummy_pos(j,3))^2);
   xi2_15(j-14)=sqrt((smooth_dummy_pos15(j+1,1)-smooth_dummy_pos15(j,1))^2+(smooth_dummy_pos15(j+1,2)-smooth_dummy_pos15(j,2))^2+(smooth_dummy_pos15(j+1,3)-smooth_dummy_pos15(j,3))^2);
end
for k=1:length(xi1_15)
    xi11_15(k)=sum(xi1_15(1:k));
end
xi1_interp15 = 0:dims(1)/2:xi11_15(end);
A_interp15 = interp1(xi11_15,A_15,xi1_interp15);
A_dash15 = gradient(A_interp15(10:end),dims(1)/2);
figure
plot(xi1_interp15(10:end),A_interp15(10:end))
disp('when the span is 15, the mean is:')
mean(abs(A_dash15))

for j=25:(length(dummy_pos)-25)
   A_25(j-24)=sqrt((dummy_pos(j,1)-smooth_dummy_pos25(j,1))^2+(dummy_pos(j,2)-smooth_dummy_pos25(j,2))^2+(dummy_pos(j,3)-smooth_dummy_pos25(j,3))^2);
   xi1_25(j-24)=sqrt((dummy_pos(j+1,1)-dummy_pos(j,1))^2+(dummy_pos(j+1,2)-dummy_pos(j,2))^2+(dummy_pos(j+1,3)-dummy_pos(j,3))^2);
   xi2_25(j-24)=sqrt((smooth_dummy_pos25(j+1,1)-smooth_dummy_pos25(j,1))^2+(smooth_dummy_pos25(j+1,2)-smooth_dummy_pos25(j,2))^2+(smooth_dummy_pos25(j+1,3)-smooth_dummy_pos25(j,3))^2);
end
for k=1:length(xi1_25)
    xi11_25(k)=sum(xi1_25(1:k));
end
xi1_interp25 = 0:dims(1)/2:xi11_25(end);
A_interp25 = interp1(xi11_25,A_25,xi1_interp25);
A_dash25 = gradient(A_interp25(10:end),dims(1)/2);
figure
plot(xi1_interp25(10:end),A_interp25(10:end))
disp('when the span is 25, the mean is:')
mean(abs(A_dash25))

for j=30:(length(dummy_pos)-30)
   A_30(j-29)=sqrt((dummy_pos(j,1)-smooth_dummy_pos30(j,1))^2+(dummy_pos(j,2)-smooth_dummy_pos30(j,2))^2+(dummy_pos(j,3)-smooth_dummy_pos30(j,3))^2);
   xi1_30(j-29)=sqrt((dummy_pos(j+1,1)-dummy_pos(j,1))^2+(dummy_pos(j+1,2)-dummy_pos(j,2))^2+(dummy_pos(j+1,3)-dummy_pos(j,3))^2);
   xi2_30(j-29)=sqrt((smooth_dummy_pos30(j+1,1)-smooth_dummy_pos30(j,1))^2+(smooth_dummy_pos30(j+1,2)-smooth_dummy_pos30(j,2))^2+(smooth_dummy_pos30(j+1,3)-smooth_dummy_pos30(j,3))^2);
end
for k=1:length(xi1_30)
    xi11_30(k)=sum(xi1_30(1:k));
end
xi1_interp30 = 0:dims(1)/2:xi11_30(end);
A_interp30 = interp1(xi11_30,A_30,xi1_interp30);
A_dash30 = gradient(A_interp30(10:end),dims(1)/2);
figure
plot(xi1_interp30(10:end),A_interp30(10:end))
disp('when the span is 30, the mean is:')
mean(abs(A_dash30))

for j=50:(length(dummy_pos)-50)
   A_50(j-49)=sqrt((dummy_pos(j,1)-smooth_dummy_pos50(j,1))^2+(dummy_pos(j,2)-smooth_dummy_pos50(j,2))^2+(dummy_pos(j,3)-smooth_dummy_pos50(j,3))^2);
   xi1_50(j-49)=sqrt((dummy_pos(j+1,1)-dummy_pos(j,1))^2+(dummy_pos(j+1,2)-dummy_pos(j,2))^2+(dummy_pos(j+1,3)-dummy_pos(j,3))^2);
   xi2_50(j-49)=sqrt((smooth_dummy_pos50(j+1,1)-smooth_dummy_pos50(j,1))^2+(smooth_dummy_pos50(j+1,2)-smooth_dummy_pos50(j,2))^2+(smooth_dummy_pos50(j+1,3)-smooth_dummy_pos50(j,3))^2);
end
for k=1:length(xi1_50)
    xi11_50(k)=sum(xi1_50(1:k));
end
xi1_interp50 = 0:dims(1)/2:xi11_50(end);
A_interp50 = interp1(xi11_50,A_50,xi1_interp50);
A_dash50 = gradient(A_interp50(10:end),dims(1)/2);
figure
plot(xi1_interp50(10:end),A_interp50(10:end))
disp('when the span is 50, the mean is:')
mean(abs(A_dash50))

for j=100:(length(dummy_pos)-100)
   A_100(j-99)=sqrt((dummy_pos(j,1)-smooth_dummy_pos100(j,1))^2+(dummy_pos(j,2)-smooth_dummy_pos100(j,2))^2+(dummy_pos(j,3)-smooth_dummy_pos100(j,3))^2);
   xi1_100(j-99)=sqrt((dummy_pos(j+1,1)-dummy_pos(j,1))^2+(dummy_pos(j+1,2)-dummy_pos(j,2))^2+(dummy_pos(j+1,3)-dummy_pos(j,3))^2);
   xi2_100(j-99)=sqrt((smooth_dummy_pos100(j+1,1)-smooth_dummy_pos100(j,1))^2+(smooth_dummy_pos100(j+1,2)-smooth_dummy_pos100(j,2))^2+(smooth_dummy_pos100(j+1,3)-smooth_dummy_pos100(j,3))^2);
end
for k=1:length(xi1_100)
    xi11_100(k)=sum(xi1_100(1:k));
end
xi1_interp100 = 0:dims(1)/2:xi11_100(end);
A_interp100 = interp1(xi11_100,A_100,xi1_interp100);
A_dash100 = gradient(A_interp100(10:end),dims(1)/2);
figure
plot(xi1_interp100(10:end),A_interp100(10:end))
disp('when the span is 100, the mean is:')
mean(abs(A_dash100))

for j=200:(length(dummy_pos)-200)
   A_200(j-199)=sqrt((dummy_pos(j,1)-smooth_dummy_pos200(j,1))^2+(dummy_pos(j,2)-smooth_dummy_pos200(j,2))^2+(dummy_pos(j,3)-smooth_dummy_pos200(j,3))^2);
   xi1_200(j-199)=sqrt((dummy_pos(j+1,1)-dummy_pos(j,1))^2+(dummy_pos(j+1,2)-dummy_pos(j,2))^2+(dummy_pos(j+1,3)-dummy_pos(j,3))^2);
   xi2_200(j-199)=sqrt((smooth_dummy_pos200(j+1,1)-smooth_dummy_pos200(j,1))^2+(smooth_dummy_pos200(j+1,2)-smooth_dummy_pos200(j,2))^2+(smooth_dummy_pos200(j+1,3)-smooth_dummy_pos200(j,3))^2);
end
for k=1:length(xi1_200)
    xi11_200(k)=sum(xi1_200(1:k));
end
xi1_interp200 = 0:dims(1)/2:xi11_200(end);
A_interp200 = interp1(xi11_200,A_200,xi1_interp200);
A_dash200 = gradient(A_interp200(10:end),dims(1)/2);
figure
plot(xi1_interp200(10:end),A_interp200(10:end))
disp('when the span is 200, the mean is:')
mean(abs(A_dash200))

clear smooth_dummy_pos20 smooth_dummy_pos5 smooth_dummy_pos10 smooth_dummy_pos15 smooth_dummy_pos25 smooth_dummy_pos30 smooth_dummy_pos50 smooth_dummy_pos100 smooth_dummy_pos200

x=[5,10,15,20,25,30,50,100,200];
y=[mean(abs(A_dash5)),mean(abs(A_dash10)),mean(abs(A_dash15)),mean(abs(A_dash20)),mean(abs(A_dash25)),mean(abs(A_dash30)),mean(abs(A_dash50)),mean(abs(A_dash100)),mean(abs(A_dash200))];
plot(x,y,'r','Linewidth',2)
xlabel('span')
ylabel('')

line_count=line_count+1;