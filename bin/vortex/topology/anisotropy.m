%script to read in the anisotropy information and plot
%if given the option print will print to .eps file rather than screen
function anisotropy(option)
if nargin==0     
  option='empty';
end
switch option
case 'print'
    disp('will not print to screen but instead to .eps files')
case 'empty'
    otherwise
    disp('incorrect option, aborting script and printing help:')
    help anisotropy
    return
end
A=load('data/anisotropy.log');
t=A(:,1) ; lpara=A(:,2) ; lperp=A(:,3);
switch option
  case 'print'
    figure('visible','off');
  otherwise
    figure('Name', 'anisotropy')      
end
  subplot(3,1,1)
    plot(t,lpara,'-r','LineWidth',2);
    set(gca,'FontSize',14)
    xlabel('t','FontSize',14)
    ylabel('lpara','FontSize',14)
  subplot(3,1,2)
    plot(t,lperp,'-b','LineWidth',2);
    set(gca,'FontSize',14)
    xlabel('t','FontSize',14)
    ylabel('lperp','FontSize',14)
  subplot(3,1,3)
    plot(t,lperp./lpara,'-k','LineWidth',2);
    set(gca,'FontSize',14)
    xlabel('t','FontSize',14)
    ylabel('lperp/lpara','FontSize',14)
switch option
  case 'print'
    disp('printing to aniso.eps')
    print('-depsc','./aniso.eps')
end

 t_Aniso=A(:,1);
save ./Anisotropy.mat t_Aniso lpara lperp
