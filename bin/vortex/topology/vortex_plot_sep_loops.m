%how many distinct vortex filaments do we have in the simulation at
%present print these loops to separate files with option knot
%e.g. vortex_loops_knotplot(1,'knot') to use with knotplot
function vortex_loops_knotplot(filenumber,option)
filename=sprintf('data/var%04d.log',filenumber);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%get the dimensions information from dims.log
dims=load('./data/dims.log');
box_size=dims(2);
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
vlength(1:4)=0.;
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
    for k=1:(length(dummy_pos)-2)
        dist=sqrt((dummy_pos(k+1,1)-dummy_pos(k,1))^2+(dummy_pos(k+1,2)-dummy_pos(k,2))^2+(dummy_pos(k+1,3)-dummy_pos(k,3))^2);
        if dist>10*dims(1)
            disp('here')
            dummy_pos(k,:)=NaN;
        end
        if k<length(dummy_pos)-2
            %vlength(l)=vlength(l)+dist;
        end
    end
    plot3(dummy_pos(:,1),dummy_pos(:,2),dummy_pos(:,3),'Color',cmap(l,:),'LineWidth',1.5)
    hold on
    
    clear dummy_pos
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
        break
    end
end

hold off
box on 
axis([-box_size/2 box_size/2 -box_size/2 box_size/2 -box_size/2 box_size/2])
daspect([1 1 1])
set(gca,'xtick',[]) ;  set(gca,'ytick',[]) ; set(gca,'ztick',[])


