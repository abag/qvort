%how many distinct vortex/magnetic filaments do we have in the simulation at 
%present
function vortex_loop_count(filenumber)
filename=sprintf('data/var%03d.log',filenumber)
fid=fopen(filename);
tline=fgetl(fid);
dummy=textscan(tline, '%f');
time=dummy{:};
tline=fgetl(fid);
dummy=textscan(tline, '%d');
number_of_particles=dummy{:};
zerocount=0;
for j=1:number_of_particles
  tline=fgetl(fid);
  dummy=textscan(tline, '%f64');
  dummy_vect=dummy{:};
  x(j)=dummy_vect(1);
  y(j)=dummy_vect(2);
  z(j)=dummy_vect(3);
  f(j)=dummy_vect(4);
  if f(j)==0
    zerocount=zerocount+1;
  end
end
f=uint16(f);
for j=1:number_of_particles
  if f(j)~=0
    next=j;
    break
  end 
end
line_count=0;
next_old=next;
counter(1:500)=0;
for l=1:500
      for i=1:number_of_particles
        line(l,i)=next;
        next=f(next);
        counter(l)=counter(l)+1;
        if next==next_old
          break   
          counter(l)
        end
        if next==0
          break
        end
      end
      line_count=line_count+1;
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
        line_count;
        break
      end
end
disp('total loop count is')
line_count
disp('brace yourself printing all the loop sizes!')
counter(1:line_count)
hist(counter(1:line_count))
xlabel('number of particles')
ylabel('frequency')
