A=load('data/KSwavenumbers.log');
wavenumber=A(:,1);
unit_vec=A(:,2:4);
subplot(1,2,1) ; plot(wavenumber,'-o')
xlabel('N')
ylabel('|K(N)|')
subplot(1,2,2) ; plot3(unit_vec(:,1),unit_vec(:,2),unit_vec(:,3),'xr','MarkerSize',10)
box on
axis equal
