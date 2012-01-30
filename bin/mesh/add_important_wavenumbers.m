function add_important_wavenumbers(intervortex)
range=ylim;
k_intervortex=2*pi/intervortex;
dumy(1)=range(1);
dumy(2)=range(2);
dumx(1:2)=k_intervortex;
hold on
plot(dumx,dumy,'k--')
dims=load('./data/dims.log');
k_box=2*pi/dims(2);
dumx(1:2)=k_box;
plot(dumx,dumy,'k-.')
k_delta=4*pi/dims(1);
dumx(1:2)=k_delta;
plot(dumx,dumy,'g--')


