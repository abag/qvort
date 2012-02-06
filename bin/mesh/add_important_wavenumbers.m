function add_important_wavenumbers(intervortex,annotate)
determine_intervortex = exist('intervortex','var');
do_annotate=exist('annotate','var');
dims=load('./data/dims.log');
if determine_intervortex==0
    A=load('./data/ts.log');
    intervortex=mean(1./sqrt(A(floor(0.5*length(A)):length(A),6)/dims(2)^3));
    disp('you have not entered intervortex spacing' )
    disp(sprintf('I have estimated it from the final half of the data, l=%f',intervortex))
end
range=ylim;
k_intervortex=2*pi/intervortex;
dumy(1)=range(1);
dumy(2)=range(2);
dumx(1:2)=k_intervortex;
hold on
plot(dumx,dumy,'k--')
if (do_annotate==1)
    text(dumx(1)+50,10^mean(log10(dumy)),'k_l','fontsize',16)
end
k_box=2*pi/dims(2);
dumx(1:2)=k_box;
plot(dumx,dumy,'k-.')
if (do_annotate==1)
    text(dumx(1)+10,10^mean(log10(dumy)),'k_D','fontsize',16)
end
k_delta=2*pi/dims(1);
dumx(1:2)=k_delta;
plot(dumx,dumy,'g--')
if (do_annotate==1)
    text(dumx(1)+100,10^mean(log10(dumy)),'k_d','fontsize',16)
end


