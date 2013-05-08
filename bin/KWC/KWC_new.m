function KWC_new(start,finish,skip)
dims=load('./data/dims.log');
cmap=colormap(jet(finish)) ;
counter=1
for i=start:skip:finish
    filename=sprintf('./data/KWC_spect%04d.log',i);
    A=load(filename);
    N=length(A); 
    A(2:N/2,2)=A(2:N/2,2)+A(end:-1:N/2+2,2);
    B(counter,:)=A(1:N/2,2);
    counter=counter+1;
    %polyfit(log(A(10:N/2-10,1)),log(A(10:N/2-10,2)),1)
    loglog(A(1:N/2,1),A(1:N/2,2),'Color',cmap(i,:))
    %loglog(A(1:N/2,1),(A(1:N/2,1).^4).*A(1:N/2,2),'Color',cmap(i,:))
    %semilogy(A(1:N/2,1),A(1:N/2,2),'Color',cmap(i,:))
    %fit=A(1:N/2,1).^-4;
    hold on
    %loglog(A(1:N/2,1),fit*1E3)
end
figure
k=A(1:N/2,1)
loglog(k,(k'.^3.66666).*mean(B),'Color','k') ; hold on 
loglog(k,(k'.^3.4).*mean(B),'Color','r')
loglog(k,(k'.^3).*mean(B),'Color','b')

end