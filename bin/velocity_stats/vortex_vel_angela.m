load velocity
%Angela's code from here on in

% Bin the (cartesian coordinate) velocity and make some plots
% Switch the binned velocity or on or off

gobin = 'on'    % 'on' will bin velocity and make a log lin plot
select = 'on'  % 'on' will find the maximum velocity and bin accordingly for each velocity
                % 'off' will use the pre-determined maximum and minimum to give even bins for vx, vy and vz
                % 'off1' uses values corresponding to bin1 selection;
                % 'off2' uses values corresponding to bin2 selection;
                % 'off3' uses values corresponding to 0.25*avdens selected;
                
switch gobin

    case 'on'
          
 switch select
     case 'on'
% start find the maximum velocity calculated        
        
maxvelx = max(ux)
minvelx = min(ux)

maxvely = max(uy)
minvely = min(uy)

maxvelz = max(uz)
minvelz = min(uz)


if abs(maxvelx)>abs(minvelx)
    mvx = maxvelx;
else
    mvx = abs(minvelx);
end

if abs(maxvely)>abs(minvely)
    mvy = maxvely;
else
    mvy = abs(minvely);
end

if abs(maxvelz)>abs(minvelz)
    mvz = maxvely;
else
    mvz = abs(minvelz);
end

nn=800;
 %  end case select on
    case 'off1'
        
        %real2d21v15 for t=5
        mvx = 23; 
        mvy = 23; 
        nn = 230;
        
        %mvx = 5;
        %mvy = 5;
        
        %nn=50;
     case 'off2'
        
         % real2d21v5out t=1
        mvx = 50;
        mvy = 50;
        
        nn=500;
        
     case 'off3'
         mvx = 8.5;
         mvy = 8.5;
         
         nn = 200;
        
 end % switch select
        
%%%%%%%%%%%%%%%%%end find the maximum velocity calculated                

velbincountx=zeros(1,2*nn+1);
velbincounty=zeros(1,2*nn+1);
velbincountz=zeros(1,2*nn+1);


for ii = 1:1:nn
vbinvaluesx(ii+nn+1) = ii*mvx/nn;
vbinvaluesy(ii+nn+1) = ii*mvy/nn;
vbinvaluesz(ii+nn+1) = ii*mvz/nn;
end

vbinvaluesx(nn+1) = 0;
vbinvaluesy(nn+1) = 0;
vbinvaluesz(nn+1) = 0;

for ii = 1:1:nn
    vbinvaluesx((nn+1)-ii)=-ii*mvx/nn;
    vbinvaluesy((nn+1)-ii)=-ii*mvy/nn;
    vbinvaluesz((nn+1)-ii)=-ii*mvz/nn;
end

lx = sum(size(ux))-1;
ly = sum(size(uy))-1;
lz = sum(size(uz))-1;


for ii = 1:1:lx
            for cc = 1:1:2*nn+1
               if cc == 1
                    if ux(1,ii)<=vbinvaluesx(cc)
                        velbincountx(cc)=velbincountx(cc)+1;
                    end
                else
                    if (ux(1,ii)>=vbinvaluesx(cc-1) && ux(1,ii)<=vbinvaluesx(cc))
                        velbincountx(cc)=velbincountx(cc)+1;
                    end
               end
            end
end
for ii = 1:1:ly
            for cc = 1:1:2*nn+1
                if cc == 1
                    if uy(1,ii)<=vbinvaluesy(cc)
                        velbincounty(cc)=velbincounty(cc)+1;
                    end
                else
                    if (uy(1,ii)>=vbinvaluesy(cc-1) && uy(1,ii)<=vbinvaluesy(cc))
                        velbincounty(cc)=velbincounty(cc)+1;
                    end
                end
            end
end
for ii = 1:1:lz
            for cc = 1:1:2*nn+1
                if cc == 1
                    if uz(1,ii)<=vbinvaluesz(cc)
                        velbincountz(cc)=velbincountz(cc)+1;
                    end
                else
                    if (uz(1,ii)>=vbinvaluesz(cc-1) && uz(1,ii)<=vbinvaluesz(cc))
                        velbincountz(cc)=velbincountz(cc)+1;
                    end
                end
            end
end

end % switch binning on or off

dvelx = vbinvaluesx(2)-vbinvaluesx(1);
dvely = vbinvaluesy(2)-vbinvaluesy(1);
dvelz = vbinvaluesz(2)-vbinvaluesz(1);

fracvelbincountx=velbincountx/(sum(velbincountx)*dvelx);
fracvelbincounty=velbincounty/(sum(velbincounty)*dvely);
fracvelbincountz=velbincountz/(sum(velbincountz)*dvelz);


% calculating the mean

sumpdfdvtotx=0;
sumpdfdvtoty=0;
sumpdfdvtotz=0;

sumpdfdvx=0;
sumpdfdvy=0;
dumpdfdvz=0;

totmeanx=0;
meanincx=0;

totmeany=0;
meanincy=0;

totmeanz=0;
meanincz=0;

for ii = 1:1:(2*nn+1)
    % for x
    sumpdfdvx=fracvelbincountx(ii)*(dvelx);
    sumpdfdvtotx=sumpdfdvtotx+sumpdfdvx;
    sumpdfdvx=0;
    % for y
    sumpdfdvy=fracvelbincounty(ii)*(dvely);
    sumpdfdvtoty=sumpdfdvtoty+sumpdfdvy;
    sumpdfdvy=0;
    % for z
    sumpdfdvz=fracvelbincountz(ii)*(dvelz);
    sumpdfdvtotz=sumpdfdvtotz+sumpdfdvz;
    sumpdfdvz=0;
    
     % calculating the mean
    meanincx = vbinvaluesx(ii)*fracvelbincountx(ii)*dvelx;
    totmeanx = totmeanx+meanincx;
    meanincx = 0;
    %
    meanincy = vbinvaluesy(ii)*fracvelbincounty(ii)*dvely;
    totmeany = totmeany+meanincy;
    meanincy = 0;
    %
    meanincz = vbinvaluesz(ii)*fracvelbincountz(ii)*dvelz;
    totmeanz = totmeanz+meanincz;
    meanincz = 0;
    
end

dvariancex = 0;
totvariancex = 0;
dvariancey = 0;
totvariancey = 0;
dvariancez = 0;
totvariancez = 0;


for ii = 1:1:(2*nn+1)
    % calculating the variance
    dvariancex = (vbinvaluesx(ii)-totmeanx)*(vbinvaluesx(ii)-totmeanx)*fracvelbincountx(ii)*dvelx;
    totvariancex = totvariancex+dvariancex;
    dvariancex = 0;
    %
   dvariancey = (vbinvaluesy(ii)-totmeany)*(vbinvaluesy(ii)-totmeany)*fracvelbincounty(ii)*dvely;
    totvariancey = totvariancey+dvariancey;
    dvariancey = 0;
    %
    dvariancez = (vbinvaluesz(ii)-totmeanz)*(vbinvaluesz(ii)-totmeanz)*fracvelbincountz(ii)*dvelz;
    totvariancez = totvariancez+dvariancez;
    dvariancez = 0;
    
end

% output statistical stuff:

totvariancex
totvariancey
totvariancez

standdevx = (totvariancex)^(0.5);
standdevy = (totvariancey)^(0.5);
standdevz = (totvariancez)^(0.5);

totmeanx
totmeany
totmeanz


sumpdfdvtotx
sumpdfdvtoty
sumpdfdvtotz

% evaluating the shape of a normal / gaussian distribution to compare:
pdfnormx=zeros(1,2*nn+1);
pdfnormy=zeros(1,2*nn+1);
pdfnormz=zeros(1,2*nn+1);

for ii = 1:1:(2*nn+1)
    pdfnormx(ii) = (1/(standdevx*((2*pi)^(0.5))))*exp(-((vbinvaluesx(ii)-totmeanx)*(vbinvaluesx(ii)-totmeanx))/(2*totvariancex));
    pdfnormy(ii) = (1/(standdevy*((2*pi)^(0.5))))*exp(-((vbinvaluesy(ii)-totmeany)*(vbinvaluesy(ii)-totmeany))/(2*totvariancey));
    pdfnormz(ii) = (1/(standdevz*((2*pi)^(0.5))))*exp(-((vbinvaluesz(ii)-totmeanz)*(vbinvaluesz(ii)-totmeanz))/(2*totvariancez));
end


switch gobin
    case 'on'

 figure(70);
 hold
 plot(vbinvaluesx,log10(velbincountx),'-d b')
 plot(vbinvaluesy,log10(velbincounty),'-v r')
 plot(vbinvaluesz,log10(velbincountz),'-* g')
 hold
 xlabel('velocity');
 ylabel('log_{10}(binned velocity count)')
 title('log lin plot psi derivatives')
 
figure(171);
hold
plot(log10(vbinvaluesx(nn+1:2*nn+1)),log10(fracvelbincountx(nn+1:2*nn+1)),'-d b')
plot(log10(vbinvaluesy(nn+1:2*nn+1)),log10(fracvelbincounty(nn+1:2*nn+1)),'-v r')
plot(log10(vbinvaluesz(nn+1:2*nn+1)),log10(fracvelbincountz(nn+1:2*nn+1)),'-* g')
hold
xlabel('log_{10}(velocity)');
ylabel('log_{10}(binned velocity count)')
title('Log Log plot of binned velocity')
% 
 figure(72);
 hold
 plot(vbinvaluesx,(velbincountx),'-d b')
 plot(vbinvaluesy,(velbincounty),'-v r')
 plot(vbinvaluesz,(velbincountz),'-* g')
 hold
 xlabel('velocity');
 ylabel('binned velocity count')
 title('Binned velocity')

figure(71);
hold
plot(vbinvaluesx/standdevx,log10(fracvelbincountx),'-o b')
plot(vbinvaluesy/standdevy,log10(fracvelbincounty),'-v r')
plot(vbinvaluesz/standdevz,log10(fracvelbincountz),'-* g')
%
plot(vbinvaluesx/standdevx,log10(pdfnormx),'. k')
plot(vbinvaluesy/standdevy,log10(pdfnormy),'-. k')
plot(vbinvaluesz/standdevz,log10(pdfnormz),'- k')
hold
xlabel('v_{i}/\sigma_{i}','fontsize',20);
ylabel('log_{10}(PDF(v_{i}))','fontsize',20);
set(gca,'FontSize',12)
axis([-10 10 -4 0.5])

figure(42);
hold
plot(vbinvaluesx,log10(fracvelbincountx),'-o b')
plot(vbinvaluesy,log10(fracvelbincounty),'-v r')
plot(vbinvaluesz,log10(fracvelbincountz),'-* g')
%
plot(vbinvaluesx,log10(pdfnormx),'. k')
plot(vbinvaluesy,log10(pdfnormy),'-. k')
plot(vbinvaluesz,log10(pdfnormz),'- k')
hold
xlabel('v_{i}','fontsize',16);
ylabel('log_{10}(PDF(v_{i}))','fontsize',16);
set(gca,'FontSize',16)
axis([-5 5 -4 0.5])

figure(12);
hold
plot(vbinvaluesy,log10(pdfnormy),'-. k','LineWidth',2)
plot(vbinvaluesx,log10(pdfnormx),': k','LineWidth',2)
plot(vbinvaluesz,log10(pdfnormz),'- k','LineWidth',2)
%
plot(vbinvaluesx,log10(fracvelbincountx),'-o b')
plot(vbinvaluesy,log10(fracvelbincounty),'-v r')
plot(vbinvaluesz,log10(fracvelbincountz),'-* g')
hold
xlabel('v_{i}','fontsize',20);
ylabel('log_{10}(PDF(v_{i}))','fontsize',20);
set(gca,'FontSize',20)
axis([-5 5 -4 0.5])

figure(700);%subplot(1,2,2)
hold
plot(log10(vbinvaluesx(nn+1:2*nn+1)),log10(fracvelbincountx(nn+1:2*nn+1)),'-o b') % b o ; m +
plot(log10(vbinvaluesy(nn+1:2*nn+1)),log10(fracvelbincounty(nn+1:2*nn+1)),'-v r') % r v ; c d
plot(log10(vbinvaluesz(nn+1:2*nn+1)),log10(fracvelbincountz(nn+1:2*nn+1)),'-* g') % g * ; y x
%plot(log10(vbinvaluesx(nn+1:2*nn+1)),-0.1422445-3.26*log10(vbinvaluesx(nn+1:2*nn+1)),'-k','Linewidth',2)
hold
xlabel('log_{10}(v_{i})','fontsize',20);
%ylabel('log_{10}(PDF(v_{i}))','fontsize',16)
set(gca,'FontSize',20);
%axis([-0.5 1.8 -4.4 -0.4])


% 
% for ii = nn+1:1:2*nn+1
%     datmat2dxr(ii-nn,1) = log10(vbinvaluesx(ii));
%     datmat2dxr(ii-nn,2) = log10(fracvelbincountx(ii));
%     datmat2dyr(ii-nn,1) = log10(vbinvaluesy(ii));
%     datmat2dyr(ii-nn,2) = log10(fracvelbincounty(ii));
% end
% 
% for ii = 1:1:nn
%     datmat2dxl(ii,1) = log10(abs(vbinvaluesx(ii)));
%     datmat2dxl(ii,2) = log10(fracvelbincountx(ii));
%     datmat2dyl(ii,1) = log10(abs(vbinvaluesy(ii)));
%     datmat2dyl(ii,2) = log10(fracvelbincounty(ii));
% end

end
