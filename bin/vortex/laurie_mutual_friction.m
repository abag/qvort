clear all ; close all ; clc
global f u u_mf v_curv v_stretch f_mf
global t_recon time
global number_of_particles
counter=0
counter2=0
for i=1:2
    vortex_load(i)
    for j=1:number_of_particles
        if f(j)==0
        else
            counter=counter+1;
            data(counter,1)=u(j);
            data(counter,2)=f_mf(j);
            if (time-t_recon(j))<0.001 && t_recon(j)>0.
                counter2=counter2+1;
                data2(counter2,1)=u(j);
                data2(counter2,2)=f_mf(j);
            end
        end
    end
end
max(v_curv)
[bandwidth,density,X,Y]=kde2d(data,2^7,[0. 0.],[max(u) max(f_mf)]);
pcolor(X,Y,density)
colormap gray
hold on
plot(data2(:,1),data2(:,2),'ro','MarkerSize',8,'MarkerFaceColor','r')
shading interp
