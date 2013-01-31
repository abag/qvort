function mesh_long_structure(filenumbers,power,fit,varargin)
% 1st input argument: input the filenumber to be analysed either single or
    % mulitple files accepted
% 2nd input argument: input the order of the structure function e.g. 2 for
    % 2nd order, or 3 for 3rd order
% 3rd input argument: either put a 1 to do a fit or a 0 to have no fit.
% 4th input argument: either type 'normal' or 'super' depending on which
    % fluid component you want to analyse
close all
optargin = size(varargin,2);
%set options based on varargin
for i=1:optargin
    switch  cell2str(varargin(i))
        case 'normal'
            normal=1;
            super=0;
        case 'super'
            super=1;
            normal=0;
        otherwise
            disp('invalid option in input arguements')
            disp('aborting code')
            return
    end
end
if fit==1
    do_fit=1;  % make this 1 if you want it to do a fit, 0 if you don't!
elseif fit==0
    do_fit=0;
else
    disp('invalid option in input argument for the fit, use either 0 or 1 in the 3rd argument')
    disp('aborting code')
    return
end
pow=power;
fluid='varargin';
load ./data/dims.log;
msize=dims(3);
if (msize==0) 
  disp('mesh size is zero exiting script')
  return
end
number_of_files=0;
for ifile=filenumbers
    number_of_files=number_of_files+1;
    filename=sprintf('data/mesh%03d.dat',ifile);
    fid=fopen(filename);
    if fid<0
        disp('mesh file does not exist, exiting script')
        return
    end
    disp(sprintf('mesh size is: %04d',msize))
    t=fread(fid,1,'float64');
    x=fread(fid,msize,'float64');
    unormx=fread(fid,msize^3,'float64');
    unormy=fread(fid,msize^3,'float64');
    unormz=fread(fid,msize^3,'float64');
    unorm_mrms=max(sqrt(unormx(:).^2+unormy(:).^2+unormz(:).^2));
    ux=fread(fid,msize^3,'float64');
    uy=fread(fid,msize^3,'float64');
    uz=fread(fid,msize^3,'float64');
    u_mrms=max(sqrt(ux(:).^2+uy(:).^2+uz(:).^2));
    unormx=reshape(unormx,msize,msize,msize);
    unormy=reshape(unormy,msize,msize,msize);
    unormz=reshape(unormz,msize,msize,msize);
    ux=reshape(ux,msize,msize,msize);
    uy=reshape(uy,msize,msize,msize);
    uz=reshape(uz,msize,msize,msize);
    if super==1
        vel_x=ux;
        vel_y=uy;
        vel_z=uz;
    end
    if normal==1
        vel_x=unormx;
        vel_y=unormy;
        vel_z=unormz;
    end
    S2(1:(msize/2))=0;
    for k=1:msize
        for j=1:msize
            for i=1:(msize/2)
                for r=1:(msize/2)
                    S2(r)=S2(r)+(abs((vel_x(k,j,i+r)-vel_x(k,j,i)))^pow);
                end
            end
        end
    end
    for k=1:msize
        for i=1:msize
            for j=1:(msize/2)
                for r=1:(msize/2)
                    S2(r)=S2(r)+(abs((vel_y(k,j+r,i)-vel_y(k,j,i)))^pow);
                end
            end
        end
    end
    for i=1:msize
        for j=1:msize
            for k=1:(msize/2)
                for r=1:(msize/2)
                    S2(r)=S2(r)+(abs((vel_z(k+r,j,i)-vel_z(k,j,i)))^pow);
                end
            end
        end
    end
end
S2=S2./(number_of_files*3*msize*msize*msize/2);
save struture.mat S2
figure('Name',strcat('2nd order long struct func., fluid:',fluid));
loglog((1:msize/2),(S2(1:msize/2)),'LineWidth',2)
xlabel('r','FontSize',14)
ylabel('struct. func.','FontSize',14)
set(gca,'Fontsize',14)
cfit = polyfit(log(1:round(msize/3)),log(S2(1:round(msize/3))),1);
if do_fit == 1
    disp(sprintf('creating a fit'))
    cfit
    hold on
    loglog((1:msize/3),exp(cfit(2)).*(1:msize/3).^(cfit(1)),'--g','LineWidth',2);
end
hold on
loglog((1:msize/2),exp(cfit(2)).*(1:msize/2).^(2/3),'--r','LineWidth',2);
hold off