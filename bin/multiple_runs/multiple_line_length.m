function multiple_line_length(nruns)
if nargin<1
 nruns=1
end
thin=4;
for i=1:nruns
    filename=sprintf('./run%d/data/ts.log',i);
    filename2=sprintf('./run%d/data/dims.log',i);
    A=load(filename);
    B=load(filename2);
    meaner(i)=mean(sqrt(A(floor(0.2*length(A)):length(A),6)/B(2)^3));
    stder(i)=std(sqrt(A(floor(0.2*length(A)):length(A),6)/B(2)^3));  
    C=A(1:thin:length(A),:);
    if (i==1)
        plot(A(:,2),A(:,6)/B(2)^3,'k-',C(:,2),C(:,6)/B(2)^3,'ko')
    elseif(i==2)
        plot(A(:,2),A(:,6)/B(2)^3,'b-',C(:,2),C(:,6)/B(2)^3,'bs')
    elseif(i==3)
        plot(A(:,2),A(:,6)/B(2)^3,'r-',C(:,2),C(:,6)/B(2)^3,'r+')
    elseif(i==4)
        plot(A(:,2),A(:,6)/B(2)^3,'g-',C(:,2),C(:,6)/B(2)^3,'g*')
    elseif(i==5)
        plot(A(:,2),A(:,6)/B(2)^3,'y-',C(:,2),C(:,6)/B(2)^3,'y^')
    elseif(i==6)
        plot(A(:,2),A(:,6)/B(2)^3,'m-',C(:,2),C(:,6)/B(2)^3,'md')
    elseif(i==7)
        plot(A(:,2),A(:,6)/B(2)^3,'c-',C(:,2),C(:,6)/B(2)^3,'cp')
    end
    hold on
end
meaner'
stder'
return
hold off
xlabel('t','FontSize',16)
ylabel('L','FontSize',16)
set(gca,'FontSize',16)
figure
for i=1:nruns
    filename=sprintf('./run%d/data/ts.log',i);
    A=load(filename);
    C=A(1:thin:length(A),:);
    for j=1:length(A)-1
        rrate(j)=(A(j+1,4)-A(j,4))/((A(j+1,2)-A(j,2))*A(j+1,6));
    end
    for j=1:length(C)-1
        rrate2(j)=(C(j+1,4)-C(j,4))/((C(j+1,2)-C(j,2))*C(j+1,6));
    end
    if (i==1)
        semilogy(C(1:length(C)-1,2),rrate2,'k-o')
    elseif(i==2)
        semilogy(C(1:length(C)-1,2),rrate2,'b-s')
    elseif(i==3)
        semilogy(C(1:length(C)-1,2),rrate2,'r-+')
    elseif(i==4)
        semilogy(C(1:length(C)-1,2),rrate2,'g-*')
    elseif(i==5)
        semilogy(C(1:length(C)-1,2),rrate2,'y-^')
    elseif(i==6)
        semilogy(C(1:length(C)-1,2),rrate2,'m-d')
    elseif(i==7)
        semilogy(C(1:length(C)-1,2),rrate2,'c-p')
    end
    hold on
    clear A C rrate2 rrate
end
xlabel('t','FontSize',16)
ylabel('R','FontSize',16)
set(gca,'FontSize',16)
return
figure
for i=1:nruns
    filename=sprintf('./run%d/data/anisotropy.log',i);
    A=load(filename);
    index = find( A(:,2) > 1.);
    if (i==1)
        plot(A(:,1),A(:,2),'k-o')
    elseif(i==2)
        plot(A(:,1),A(:,2),'b-s')
    elseif(i==3)
        plot(A(:,1),A(:,2),'r-+')
    elseif(i==4)
        plot(A(:,1),A(:,2),'g-*')
    elseif(i==5)
        plot(A(:,1),A(:,2),'y-^')
    elseif(i==6)
        plot(A(:,1),A(:,2),'m-d')
    elseif(i==7)
        plot(A(:,1),A(:,2),'c-p')
    end
    hold on
    clear index A
end
hold off
xlabel('t','FontSize',16)
ylabel('l','FontSize',16)
set(gca,'FontSize',16)
