function multiple_line_length(nruns)
if nargin<1
 nruns=1
end
for i=1:nruns
    filename=sprintf('./run%d/data/ts.log',i);
    filename2=sprintf('./run%d/data/dims.log',i);
    A=load(filename);
    B=load(filename2);
    meaner(i)=mean(A(floor(0.2*length(A)):length(A),6)/B(2)^3);
    stder(i)=std(A(floor(0.2*length(A)):length(A),6)/B(2)^3);  
    if (i==1)
        plot(A(:,2),A(:,6)/B(2)^3,'k-o')
    elseif(i==2)
        plot(A(:,2),A(:,6)/B(2)^3,'b-s')
    elseif(i==3)
        plot(A(:,2),A(:,6)/B(2)^3,'r-+')
    elseif(i==4)
        plot(A(:,2),A(:,6)/B(2)^3,'g-*')
    elseif(i==5)
        plot(A(:,2),A(:,6)/B(2)^3,'y-^')
    elseif(i==6)
        plot(A(:,2),A(:,6)/B(2)^3,'m-d')
    elseif(i==7)
        plot(A(:,2),A(:,6)/B(2)^3,'c-p')
    end
    hold on
end
%meaner'
%stder'
%return
hold off
xlabel('t','FontSize',16)
ylabel('L','FontSize',16)
set(gca,'FontSize',16)
figure
for i=1:nruns
    filename=sprintf('./run%d/data/ts.log',i);
    A=load(filename);
    for j=1:length(A)-1
        rrate(j)=(A(j+1,4)-A(j,4))/((A(j+1,2)-A(j,2))*A(j+1,6));
    end
    if (i==1)
        semilogy(A(1:length(A)-1,2),rrate,'k-o')
    elseif(i==2)
        semilogy(A(1:length(A)-1,2),rrate,'b-s')
    elseif(i==3)
        semilogy(A(1:length(A)-1,2),rrate,'r-+')
    elseif(i==4)
        semilogy(A(1:length(A)-1,2),rrate,'g-*')
    elseif(i==5)
        semilogy(A(1:length(A)-1,2),rrate,'y-^')
    elseif(i==6)
        semilogy(A(1:length(A)-1,2),rrate,'m-d')
    elseif(i==7)
        semilogy(A(1:length(A)-1,2),rrate,'c-p')
    end
    hold on
    clear A rrate
end
xlabel('t','FontSize',16)
ylabel('R','FontSize',16)
set(gca,'FontSize',16)
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
