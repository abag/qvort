function density_box_measure(filenumbers,ndivide)
%load in time series data to get number of files
A=load('./data/ts.log');
%get the dimensions information from dims.log
dims=load('./data/dims.log');
counter=1
for i=filenumbers
    i
    filename=sprintf('data/var%04d.log',i);
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
        fclose(fid);
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
        end
        f=uint16(f);
        fclose(fid);
    end
    for l=0:ndivide-1
        for m=0:ndivide-1
            for n=0:ndivide-1
                 total_length(counter)=0. ;
                min_boundsx=-dims(2)/2.+l*dims(2)/ndivide;
                max_boundsx=-dims(2)/2.+(l+1)*dims(2)/ndivide;
                min_boundsy=-dims(2)/2.+m*dims(2)/ndivide;
                max_boundsy=-dims(2)/2.+(m+1)*dims(2)/ndivide;
                min_boundsz=-dims(2)/2.+n*dims(2)/ndivide;
                max_boundsz=-dims(2)/2.+(n+1)*dims(2)/ndivide;
                for j=1:number_of_particles
                    if round(f(j))==0
                    else
                        dummy_x(1,1)=x(j);
                        dummy_x(2,1)=x(round(f(j)));
                        dummy_x(1,2)=y(j);
                        dummy_x(2,2)=y(round(f(j)));
                        dummy_x(1,3)=z(j);
                        dummy_x(2,3)=z(round(f(j)));
                        dummy_x;
                        if dummy_x(1,1)>min_boundsx && dummy_x(1,1)<max_boundsx && dummy_x(1,2)>min_boundsy && dummy_x(1,2)<max_boundsy && dummy_x(1,3)>min_boundsz && dummy_x(1,3)<max_boundsz
                            dist=sqrt((dummy_x(1,1)-dummy_x(2,1))^2+(dummy_x(1,2)-dummy_x(2,2))^2+(dummy_x(1,3)-dummy_x(2,3))^2);
                            if (dist<dims(2)/4.)
                                total_length(counter)=total_length(counter)+dist;
                                %total_length
                            end
                        end
                    end
                end
                counter=counter+1;
            end
        end
    end
end
ksdensity(total_length/(dims(2)/ndivide)^3)
m_vol=mean(total_length/(dims(2)/ndivide)^3);
sd_vol=std(total_length/(dims(2)/ndivide)^3);
dummy_x=linspace(min(total_length/(dims(2)/ndivide)^3),max(total_length/(dims(2)/ndivide)^3),100);
hold on
plot(dummy_x,normpdf(dummy_x,m_vol,sd_vol),'r');
