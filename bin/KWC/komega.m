function komega(filenumber)
filename=sprintf('./data/komega%04d.dat',filenumber);
fid=fopen(filename);
 A=fread(fid,'float64');
 size(A);
 fclose(fid)
 A=reshape(A,2,327,256);
 B=squeeze(A(1,1:300,:))+sqrt(-1)*squeeze(A(2,1:300,:));
 B=B';
 %figure
 %pcolor(abs(B)) ; shading interp
 C=fft2(B);
 n1=256;
 n2=300;
 C(1,:)=[];
 k1= [-n1/2:-1,1:n1/2-1];
 k2= (2/n2)*[-n2/2:-1,0:n2/2-1]*2*pi/(14*1E-2);
 %figure
 D=fftshift(C');
 size(D)
 size(k1)
 for i=1:size(D,1)
     k11=k1.^4;
     D(i,:)=D(i,:).*((k11));
 end
 pcolor(k1,k2,log(abs(D))) ; shading interp
 %pcolor(k1,k2,(abs(fftshift(C')))) ; shading interp
 hold on
 CC=2*pi %log(2)-0.57
 plot(k1,(1E-3)*((k1.^2)/(4*pi)).*(log(1./(abs(k1)*8.244023E-9))+CC),'r')
 %plot(k1,(9.97E-4)*((k1.^2)/(4*pi)).*(log(1./(abs(k1)*1E-10))+CC),'b')
 %plot(k1,(9.97E-4)*((k1.^2)/(4*pi)).*(log(1./(abs(k1)*1E-10))+CC),'b')
% plot(k1,(1.3E-3)*((k1.^2)),'b')
