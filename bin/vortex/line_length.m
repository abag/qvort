%read in the ts file and plot the line length in both normal and logscale
%takes an input that is the number of sections to divide the plot into and give
%a growth/decay rate
function line_length(sections)
if nargin<1
  sections=1;
end
A=load('./data/ts.log');
dims=load('./data/dims.log');
t=A(:,2) ; l=A(:,6)/dims(2)^3 ; 
scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4) scrsz(3)/1.2 scrsz(4)/2],'PaperPosition',[0.25 2.5 28.0 12.0],'color','w','visible','on','Name', 'line_length')
  subplot(1,2,1)
    plot(t,l,'-r','LineWidth',2);
    set(gca,'FontSize',14)
    xlabel('t','FontSize',14)
    ylabel('l','FontSize',14)
  subplot(1,2,2)
    plot(t,log(l),'-b','LineWidth',2);
    hold on
    for i=1:sections
      start=floor((i-1)/sections*length(t)+1);
      finish=ceil(i/sections*length(t));
      p=polyfit((t(start:finish)),log(l(start:finish)),1);
      dummy_ll=p(1)*(t(start:finish))+p(2);
      plot((t(start:finish)),dummy_ll,'-.r','LineWidth',2);
      text((t(floor(start/2+finish/2))),dummy_ll(floor(length(dummy_ll)/2))-0.05*(max(log(l))-min(log(l))),num2str(p(1)),'FontSize',14)
    end
    hold off
    xlabel('t','FontSize',14)
    ylabel('log l','FontSize',14)
    set(gca,'FontSize',14)
figure('Name', 'line density only')
    plot(t,l,'-k','LineWidth',2);
    set(gca,'FontSize',16)
    xlabel('t','FontSize',16)
    ylabel('L','FontSize',16)
