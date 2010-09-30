function varargout = showvol(varargin)
%  Function Showvol, for showing an phong-shaded iso-surface of volume data.
%     showvol (Volume, isovalue) 
%   or
%     showvol (Volume)
%  
% example:
%  % Load Vessel Data
%   load('CT_aneurysm');
%  % Show iso surface
%   showvol(A);
%
% Function is written by D.Kroon University of Twente (November 2009)


% Last Modified by GUIDE v2.5 30-Nov-2009 15:29:49

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @showvol_OpeningFcn, ...
                   'gui_OutputFcn',  @showvol_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before showvol is made visible.
function showvol_OpeningFcn(hObject, eventdata, handles, varargin)

% Choose default command line output for showvol
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes showvol wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% Check if input volume is specified
if (isempty(varargin)), error('showgl:inputs', 'no input image'); end

% Get the input voxel volume
data=struct;
data.volume=varargin{1};

% Calculate min and max value to scale the intensities
data.vmax=double(max(data.volume(:)));
data.vmin=double(min(data.volume(:)));

% Get or Calculate the staritng iso value
if(length(varargin)>1), data.iso=double(varargin{2});
else data.iso=(3/4)*(data.vmax-data.vmin)-data.vmin; end

% Default input values
data.scaling=0.5;
data.posx=0.5; 
data.posy=0.5;  
data.posz=0.5; 
data.scales=[1 1 1];
data.color=[1 0.75 0.65];

data.sizes=size(data.volume);
data.handles=handles;
data.isorender=false;
data.slicerender=false;
data.h1=[]; 
data.h2=[]; 
data.h3=[];

% Store all data in the Gui
setMyData(data);

% Make downsampled versions for fast rendering
make_data_small();

% Create and show a log histogram of the data
make_histogram();
show_histogram();

% Show an preview of the iso surface
show_preview()

% This function creates an histogram image of the voxel data
function make_histogram()
data=getMyData();
histdata=histc(data.volume(:),linspace(data.vmin,data.vmax,350));
histdata=log(histdata+1);
histdata=histdata./max(histdata(:));
histimage=zeros(100,350);
for i=1:350, histimage(round(1:double(histdata(i).*100)),i)=1; end
data.histimage=histimage(end:-1:1,:);
setMyData(data);

% This function shows the histogram image, with a part green and part
% red, based on the current iso value
function show_histogram()
data=getMyData();
I=zeros([size(data.histimage) 3]);
j=round(size(data.histimage,2)*(data.iso-data.vmin)/(data.vmax-data.vmin));
j(j<1)=1;
I(:,1:j,1)=data.histimage(:,1:j);
I(:,j:end,2)=data.histimage(:,j:end);
imshow(I,'Parent',data.handles.axes_histogram);
set(data.handles.slider_iso,'value',(data.iso-data.vmin)/(data.vmax-data.vmin))
set(data.handles.edit3,'String',num2str(data.iso));

% This function makes downsampled volumes, for fast iso surface rendering. 
function make_data_small()
data=getMyData();
    data.small1=imresize3d(data.volume,0.25,[],'linear');
    data.small2=imresize3d(data.volume,0.5,[],'linear');
setMyData(data);

% Make the vertex, face and normals list of the iso-surface
function make_isosurface()
data=getMyData();
if(data.scaling==0.25)
    data.FV = isosurface(data.small1, data.iso);
    data.N  = isonormals(smooth3(data.small1),data.FV.vertices);
elseif(data.scaling==0.5)
    data.FV = isosurface(data.small2, data.iso);
    data.N  = isonormals(smooth3(data.small2),data.FV.vertices);
else
    data.FV = isosurface(data.volume, data.iso);
    data.N  = isonormals(smooth3(data.volume),data.FV.vertices);
end
data.FV.vertices=data.FV.vertices*(1/data.scaling);
data.FVmean=mean(data.FV.vertices,1);
setMyData(data);

% This function is a small wraper around the render_iso function,
% which uses some basic Matlab functions to render within 10th of a
% second a iso-surface to an image.
function show_preview()
data=getMyData();
J=render_iso(data.volume,1,data.iso);
J=imrotate90(J);
imshow(J,'Parent',data.handles.axes_preview);
setMyData(data);


% This function displays the iso-surface polygon data
function show_isosurface()
data=getMyData();
cla(data.handles.axes3d);
data.isorender=true;
data.handle_patch= patch(data.FV, 'Facecolor', data.color, 'EdgeColor', 'none','VertexNormals',data.N);
material shiny
data.handel_light=camlight('headlight');
axis equal; daspect(data.scales); axis off; 
setMyData(data);


% --- Outputs from this function are returned to the command line.
function varargout = showvol_OutputFcn(hObject, eventdata, handles) 
% Get default command line output from handles structure
varargout{1} = handles.output;


function setMyData(data)
% Store data struct in figure
setappdata(gcf,'data3d',data);

function data=getMyData()
% Get data struct stored in figure
data=getappdata(gcf,'data3d');

function A=imresize3d(V,scale,tsize,ntype,npad)
% This function resizes a 3D image volume to new dimensions
% Vnew = imresize3d(V,scale,nsize,ntype,npad);
%
% inputs,
%   V: The input image volume
%   scale: scaling factor, when used set tsize to [];
%   nsize: new dimensions, when used set scale to [];
%   ntype: Type of interpolation ('nearest', 'linear', or 'cubic')
%   npad: Boundary condition ('replicate', 'symmetric', 'circular', 'fill', or 'bound')  
%
% outputs,
%   Vnew: The resized image volume
%
% example,
%   load('mri','D'); D=squeeze(D);
%   Dnew = imresize3d(D,[],[80 80 40],'nearest','bound');
%
% This function is written by D.Kroon University of Twente (July 2008)

% Check the inputs
if(exist('ntype', 'var') == 0), ntype='nearest'; end
if(exist('npad', 'var') == 0), npad='bound'; end
if(exist('scale', 'var')&&~isempty(scale)), tsize=round(size(V)*scale); end
if(exist('tsize', 'var')&&~isempty(tsize)),  scale=(tsize./size(V)); end

% Make transformation structure   
T = makehgtform('scale',scale);
tform = maketform('affine', T);

% Specify resampler
R = makeresampler(ntype, npad);

% Resize the image volueme
A = tformarray(V, tform, R, [1 2 3], [1 2 3], tsize, [], 0);


% This function is used to show cross sections of the volume
function ShowSlices
showslicex();
showslicey();
showslicez();
data=getMyData();
if(~data.slicerender), 
    data.slicerender=true; 
    view(3); axis square; axis off; 
end

setMyData(data)

function showslicex
data=getMyData();
if(ishandle(data.h1)), delete(data.h1); end
posx=max(min(round(data.sizes(1)*data.posx),data.sizes(1)),1);

if(ndims(data.volume)==3)
    slicexg = uint8(255*(double(squeeze(data.volume(posx,:,:)))-data.vmin)/(data.vmax-data.vmin))';
    slicex=ind2rgb(slicexg,gray(256));
else
    slicexs = uint8(255*(double(squeeze(data.volume(posx,:,:,:)))-data.vmin)/(data.vmax-data.vmin)); 
    slicex=uint8(zeros([size(slicexs,2) size(slicexs,1) 3]));
    for i=1:3, slicex(:,:,i)=slicexs(:,:,i)'; end
end
slicex_y=[posx posx;posx posx];
slicex_x=[0 (data.sizes(2)-1);0 (data.sizes(2)-1)];
slicex_z=[0 0;(data.sizes(3)-1) (data.sizes(3)-1)];

data.h1=surface(slicex_x,slicex_y,slicex_z, slicex,'FaceColor','texturemap', 'EdgeColor',[1 0 0], 'CDataMapping','direct','FaceAlpha',1);
setMyData(data)

function showslicey()
data=getMyData();
if(ishandle(data.h2)),delete(data.h2);  end
posy=max(min(round(data.sizes(2)*data.posy),data.sizes(2)),1); 

if(ndims(data.volume)==3)
    sliceyg = uint8(255*(double(squeeze(data.volume(:,posy,:)))-data.vmin)/(data.vmax-data.vmin))';
    slicey=ind2rgb(sliceyg,gray(256));
else
    sliceys = uint8(255*(double(squeeze(data.volume(:,posy,:,:)))-data.vmin)/(data.vmax-data.vmin)); 
    slicey=uint8(zeros([size(sliceys,2) size(sliceys,1) 3]));
    for i=1:3, slicey(:,:,i)=sliceys(:,:,i)';end
end

slicey_y=[0 (data.sizes(1)-1);0 (data.sizes(1)-1)];
slicey_x=[posy posy;posy posy]; 
slicey_z=[0 0;(data.sizes(3)-1) (data.sizes(3)-1)];

data.h2=surface(slicey_x,slicey_y,slicey_z, slicey,'FaceColor','texturemap', 'EdgeColor',[1 0 0], 'CDataMapping','direct','FaceAlpha',1);
setMyData(data)

function showslicez()
data=getMyData();
if(ishandle(data.h3)), delete(data.h3); end
posz=max(min(round(data.sizes(3)*data.posz),data.sizes(3)),1); 

if(ndims(data.volume)==3)
    slicezg = uint8(255*(double(squeeze(data.volume(:,:,posz)))-data.vmin)/(data.vmax-data.vmin));
    slicez=ind2rgb(slicezg,gray(256));
else
    slicezs = uint8(255*(double(squeeze(data.volume(:,:,posz,:)))-data.vmin)/(data.vmax-data.vmin)); 
    slicez=uint8(zeros([size(slicezs,2) size(slicezs,1) 3]));
end
slicez_x=[0 (data.sizes(1)-1);0 (data.sizes(1)-1)]; slicez_y=[0 0;(data.sizes(2)-1) (data.sizes(2)-1)]; slicez_z=[posz posz;posz posz];
data.h3=surface(slicez_x,slicez_y,slicez_z, slicez,'FaceColor','texturemap', 'EdgeColor',[1 0 0], 'CDataMapping','direct','FaceAlpha',1);
setMyData(data)

% The mouse movement function is used to update the position of the camera light.
function figure1_WindowButtonMotionFcn(hObject, eventdata, handles)
data=getMyData();
if(data.isorender)
    camlight(data.handel_light,'headlight');
end

% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function figure1_WindowButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function figure1_WindowButtonUpFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_file_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_render_Callback(hObject, eventdata, handles)
% hObject    handle to menu_render (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_iso_surface_Callback(hObject, eventdata, handles)
% hObject    handle to menu_iso_surface (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function menu_slice_view_Callback(hObject, eventdata, handles)
% hObject    handle to menu_slice_view (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on slider movement.
function slider_iso_Callback(hObject, eventdata, handles)
% hObject    handle to slider_iso (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
data=getMyData();
data.iso=get(data.handles.slider_iso,'value')*(data.vmax-data.vmin)-data.vmin;
setMyData(data);
show_preview();
show_histogram();

% --- Executes during object creation, after setting all properties.
function slider_iso_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_iso (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% Callback from iso-value edit box
function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double
data=getMyData();
data.iso=str2double(get(data.handles.edit3,'String'));
setMyData(data);
show_preview();
show_histogram();




% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% Make and render an iso surface after user pushed the button
function pushbutton_render_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_render (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data=getMyData();
axes(data.handles.axes3d); 
cla(data.handles.axes3d,'reset'); 
rotate3d('on');
data.slicerender=false;
setMyData(data);
make_isosurface();
show_isosurface();


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over slider_iso.
function slider_iso_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to slider_iso (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


function J=render_iso(K,dim,iso,Options)
% This function renders a iso surface of an image-volume to an 2D image
% 
% J = render_iso(K,dim,iso,Options)
%
% Example,
%   load MRI
%   D = squeeze(D); D=D(:,:,end:-1:1);
%   Ds= smooth3(D);
%   figure, render_iso(D,3,30,struct('S',Ds,'color',[1 0 0]));
%
%
% Function is written by D.Kroon University of Twente (November 2009)

% Check inputs
if(nargin<1), error('iso_surface:inputs','no input volume '); end
if(nargin<2), dim=1; end
if(nargin<3), iso=max(K(:))/3; end

% Process input options
defaultoptions=struct('V',[0 0 1],'L',[0.67 0.33 0.67],'ka',0.3,'kd',0.5,'ks',0.5,'n',9,'color',[0 0 1],'S',[]);
if(~exist('Options','var')),
    Options=defaultoptions;
else
    tags = fieldnames(defaultoptions);
    for i=1:length(tags)
        if(~isfield(Options,tags{i})),  Options.(tags{i})=defaultoptions.(tags{i}); end
    end
    if(length(tags)~=length(fieldnames(Options))),
        warning('render_iso:unknownoption','unknown options found');
    end
end

% Get the iso-surface location
[surfacefound,I]=max(K>iso,[],dim);
if(isempty(Options.S)), A=K; else A=Options.S; end
% Calculate image gradients on the iso-surface location with central
% differences
if(dim==1)
    [x,y]=ndgrid(1:size(A,2),1:size(A,3)); x=x(:); y=y(:);
    z=circshift(I,[0 1 0]); z=z(:);  Ixp=sub2ind(size(A),z,x,y);
    z=circshift(I,[0 -1 0]); z=z(:); Ixm=sub2ind(size(A),z,x,y);
    z=circshift(I,[0 0 1]); z=z(:);  Iyp=sub2ind(size(A),z,x,y);
    z=circshift(I,[0 0 -1]); z=z(:); Iym=sub2ind(size(A),z,x,y);
    z=min(I+1,size(A,1)); z=z(:);    Izp=sub2ind(size(A),z,x,y); 
    z=max(I-1,1); z=z(:);            Izm=sub2ind(size(A),z,x,y);
    Axp=circshift(reshape(A(Ixp),size(surfacefound)),[0 -1 0]);
    Axm=circshift(reshape(A(Ixm),size(surfacefound)),[0 1 0]);
    Ayp=circshift(reshape(A(Iyp),size(surfacefound)),[0 0 -1]);
    Aym=circshift(reshape(A(Iym),size(surfacefound)),[0 0 1] );
elseif(dim==2)
    [x,y]=ndgrid(1:size(A,1),1:size(A,3)); x=x(:); y=y(:);
    z=circshift(I,[1 0 0]); z=z(:);  Ixp=sub2ind(size(A),x,z,y);
    z=circshift(I,[-1 0 0]); z=z(:); Ixm=sub2ind(size(A),x,z,y);
    z=circshift(I,[0 0 1]); z=z(:);  Iyp=sub2ind(size(A),x,z,y);
    z=circshift(I,[0 0 -1]); z=z(:); Iym=sub2ind(size(A),x,z,y);
    z=min(I+1,size(A,1)); z=z(:);    Izp=sub2ind(size(A),x,z,y); 
    z=max(I-1,1); z=z(:);            Izm=sub2ind(size(A),x,z,y);
    Axp=circshift(reshape(A(Ixp),size(surfacefound)),[-1 0 0]);
    Axm=circshift(reshape(A(Ixm),size(surfacefound)),[1 0 0]);
    Ayp=circshift(reshape(A(Iyp),size(surfacefound)),[0 0 -1]);
    Aym=circshift(reshape(A(Iym),size(surfacefound)),[0 0 1] );
else
    [x,y]=ndgrid(1:size(A,1),1:size(A,2)); x=x(:); y=y(:);
    z=circshift(I,[1 0 0]); z=z(:);  Ixp=sub2ind(size(A),x,y,z);
    z=circshift(I,[-1 0 0]); z=z(:); Ixm=sub2ind(size(A),x,y,z);
    z=circshift(I,[0 1 0]); z=z(:);  Iyp=sub2ind(size(A),x,y,z);
    z=circshift(I,[0 -1 0]); z=z(:); Iym=sub2ind(size(A),x,y,z);
    z=min(I+1,size(A,1)); z=z(:);    Izp=sub2ind(size(A),x,y,z); 
    z=max(I-1,1); z=z(:);            Izm=sub2ind(size(A),x,y,z);
    Axp=circshift(reshape(A(Ixp),size(surfacefound)),[-1 0 0]);
    Axm=circshift(reshape(A(Ixm),size(surfacefound)),[1 0 0]);
    Ayp=circshift(reshape(A(Iyp),size(surfacefound)),[0 -1 0]);
    Aym=circshift(reshape(A(Iym),size(surfacefound)),[0 1 0] );
end
Azp=reshape(A(Izp),size(surfacefound));
Azm=reshape(A(Izm),size(surfacefound));

surfacefound=squeeze(surfacefound);
% Normalize the iso-surface imagegradient to get iso-normals
Fx=double(Axp-Axm); Fy=double(Ayp-Aym); Fz=double(Azp-Azm);
Lf=1./(sqrt(Fx.^2+Fy.^2+Fz.^2)+eps);
N=zeros([size(surfacefound) 3]);
N(:,:,1)=Fx.*Lf; N(:,:,2)=Fy.*Lf; N(:,:,3)=Fz.*Lf;

% Calculate the phong shading
L=permute(Options.L,[3 1 2]); L=repmat(L,[size(N,1) size(N,2) 1]);
V=permute(Options.V,[3 1 2]); V=repmat(V,[size(N,1) size(N,2) 1]);
Ia=Options.ka*surfacefound;
Id=Options.kd*sum(L.*N,3).*surfacefound;
t=max(sum((2.0*repmat(sum(N.*L,3),[1 1 3]).*N - L).*V,3),0);
Is=Options.ks*(t.^Options.n).*surfacefound;

% Make the final iso-surface image
J=zeros([size(surfacefound) 3]);
J(:,:,1)=Is+(Ia+Id)*Options.color(1);
J(:,:,2)=Is+(Ia+Id)*Options.color(2);
J(:,:,3)=Is+(Ia+Id)*Options.color(3);

% Show the image if there is no output variable
if(nargout == 0), imshow(J); end

     


% --- Executes on slider movement.
function sliderx_Callback(hObject, eventdata, handles)
% hObject    handle to sliderx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
data=getMyData();
data.posx=get(data.handles.sliderx,'value');
setMyData(data);
%if(data.slicerender)
    showslicex();
%end

% --- Executes during object creation, after setting all properties.
function sliderx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes on slider movement.
function slidery_Callback(hObject, eventdata, handles)
% hObject    handle to slidery (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
data=getMyData();
data.posy=get(data.handles.slidery,'value');
setMyData(data);
%if(data.slicerender)
    showslicey();
%end

% --- Executes during object creation, after setting all properties.
function slidery_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slidery (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function sliderz_Callback(hObject, eventdata, handles)
% hObject    handle to sliderz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
data=getMyData();
data.posz=get(data.handles.sliderz,'value');
setMyData(data);
%if(data.slicerender)
    showslicez();
%end

% --- Executes during object creation, after setting all properties.
function sliderz_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider11_Callback(hObject, eventdata, handles)
% hObject    handle to slider_iso (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_iso (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_render.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_render (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on button press in pushbutton_render_slices.
function pushbutton_render_slices_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_render_slices (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data=getMyData();
axes(data.handles.axes3d);
cla(data.handles.axes3d,'reset');
rotate3d('on');
data.isorender=false;
setMyData(data);
ShowSlices;

function J=imrotate90(I)
J=zeros(size(I,2),size(I,1),size(I,3));
for i=1:size(I,3), J(:,:,i)=I(:,end:-1:1,i)'; end

% --- Executes when selected object is changed in uipanel2.
function uipanel2_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel2 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
data=getMyData();
sel=[get(data.handles.radiobutton025,'value') get(data.handles.radiobutton050,'value') get(data.handles.radiobutton100,'value')];
if(sel(1))
    data.scaling=0.25;
elseif(sel(2))
    data.scaling=0.5;
else
    data.scaling=1;
end
setMyData(data);
