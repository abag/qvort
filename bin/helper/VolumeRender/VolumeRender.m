function varargout = VolumeRender(varargin)
% This function is a Volume Render application, using OpenGL surfaces. 
%
% color and alpha maps can be changed on the fly by dragging and creating
% new color/alpha markers with the left mouse button.
%
% VolumeRender(V,Scales,PosData,AlphaData,ColorData)
%
% V = 3D Volume
% Scales = Scaling of dimensions (see dicom), defaults [1 1 1]
%
% Histogram parameters, (specify ALL or NONE of the parameters below )
% PosData = Position (grey value) of markers (range [0 1]), defaults [0.2 0.4 0.6 0.9]
% AlphaData = Alpha (transparency) of markers (range [0 1]), defaults [0 0.5 0.35 1]
% ColorData = Color of markers (range [0 0 0; 1 1 1]), defaults [0 0 0; 1 0 0; 1 1 0; 1 1 1]
%
% Example:
% load sampledata;
% VolumeRender(A);
%
% Known bug: changes with function colormap consume large amount of memory.
%
% Author D.Kroon of University of Twente 2008-11-03

% Edit the above text to modify the response to help VolumeRender

% Last Modified by GUIDE v2.5 11-Mar-2008 13:32:24

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @VolumeRender_OpeningFcn, ...
                   'gui_OutputFcn',  @VolumeRender_OutputFcn, ...
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


% --- Executes just before VolumeRender is made visible.
function VolumeRender_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to VolumeRender (see VARARGIN)
% Choose default command line output for VolumeRender
    handles.output = hObject;
    % Update handles structure
    guidata(hObject, handles);

    % For detection of mouse up/down during mousedrag
    set(gcf, 'userdata', {0});
    set(gcf, 'WindowButtonDownFcn', 'userdata=get(gcf, ''userdata'');userdata{1}=1;set(gcf, ''userdata'', userdata)');
    set(gcf, 'WindowButtonUpFcn'  , 'userdata=get(gcf, ''userdata'');userdata{1}=0;set(gcf, ''userdata'', userdata)');

    % Set render to opengl
    set(gcf, 'Renderer', 'opengl'); 
    
    % Process the Inputs
    if(~isempty(varargin)), histo.Voxelvolume=varargin{1}; else histo.Voxelvolume=zeros(2,2,2); end
    if(length(varargin)>1), histo.Scales=varargin{2}; else histo.Scales=[1 1 1]; end
    if(length(varargin)>2), histo.value=varargin{3}; else histo.value = [0.2 0.4 0.6 0.9]; end
    if(length(varargin)>3), histo.alpha=varargin{4}; else histo.alpha = [0 0.5 0.35 1]; end
    if(length(varargin)>4), histo.colors=varargin{5}; else histo.colors = [0 0 0; 1 0 0; 1 1 0; 1 1 1]; end

    % Store range of voxel data
    rangc=getrangefromclass(histo.Voxelvolume);
    histo.minx=rangc(1); histo.maxx=rangc(2);
    
    % Save Handles of histogram window
    histo.histogramwindowhandle=gcf;
    histo.histogramaxeshandle=handles.axes_histogram;

    % Create Window for rendering
    figure,
    histo.mainwindowhandle=gcf;
    set(gcf, 'Renderer', 'opengl'); 
    set(gcf, 'BackingStore', 'off');
    
    % Save handles of render window
    histo.mainwindowaxeshandle=gca;
    
    figure(histo.histogramwindowhandle);
    % Create and show the histogram 
    histo=createHistogram(histo);
    % Create the markers and curve
    histo=drawPoints(histo,handles);
    % Create the color and alpha maps from the markes
    histo=createAlphaColorMap(histo);
    % Render the 3D volume
    histo=volumerender(histo);

    % Save the data in the figure
setMyData(histo);


% --- Outputs from this function are returned to the command line.
function varargout = VolumeRender_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


 
function histo=createAlphaColorMap(histo)
% This function creates a Matlab colormap and alphamap from the markers
    histo.ColorMap=zeros(255,3); histo.AlphaMap=zeros(255,1);
    % Loop through all 256 color/alpha indexes
    for j=0:255 
        i=j/255;
        if (i<histo.value(1)), alpha=0; color=histo.colors(1,:);
        elseif(i>histo.value(end)), alpha=0; color=histo.colors(end,:);
        elseif(i==histo.value(1)), alpha=histo.value(1); color=histo.colors(1,:);
        elseif(i==histo.value(end)), alpha=histo.value(end); color=histo.colors(end,:);
        else
            % Linear interpolate the color and alpha between markers
            index_down=find(histo.value<=i); index_down=index_down(end);
            index_up=find(histo.value>i); index_up=index_up(1);
            perc=(i-histo.value(index_down)) / (histo.value(index_up) - histo.value(index_down));
            color=(1-perc)*histo.colors(index_down,:)+perc*histo.colors(index_up,:);
            alpha=(1-perc)*histo.alpha(index_down)+perc*histo.alpha(index_up);
        end
        histo.AlphaMap(j+1)=alpha;
        histo.ColorMap(j+1,:)=color;
    end

    
function histo=createHistogram(histo)
% This function creates and show the (log) histogram of the data    
    % Focus on histogram axes
    axes(histo.histogramaxeshandle);
    % Get histogram
    [histo.countsy, histo.countsx]=imhist(histo.Voxelvolume(:));
    % Log the histogram data
    histo.countsy=log(histo.countsy+100); histo.countsy=histo.countsy-min(histo.countsy);
    % Display the histogram
    stem(histo.countsx,histo.countsy,'Marker', 'none'); hold on; 
    % Set the axis of the histogram axes
    histo.maxy=max(histo.countsy);
    axis([histo.minx histo.maxx 0 histo.maxy]);

function histo=volumerender(histo)
    axes(histo.mainwindowaxeshandle);
    % Initialize the figure
    hold on; axis equal; axis xy; axis off; 
    
    % Set the data dimension scaling
    set(histo.mainwindowaxeshandle,'DataAspectRatio',histo.Scales);
    
    % Set the colors and alphamap 
    colormap(histo.mainwindowaxeshandle,histo.ColorMap); 
    alphamap(histo.mainwindowhandle,histo.AlphaMap);

    sizes = size(histo.Voxelvolume);

    % Making the X slices (surfaces)
    for posx = 1:sizes(1)
        % Position of surface
        slicex_x = [posx posx;posx posx]; slicex_y = [0 (sizes(2)-1);0 (sizes(2)-1)]; slicex_z = [0 0;(sizes(3)-1) (sizes(3)-1)];

        % Texture of surface 
        slicex = squeeze(histo.Voxelvolume(posx,:,:)); slicex = im2uint8(slicex)';

        % Transparance of surface
        maskx = slicex;

        % Display Surface
        surface(slicex_x,slicex_y,slicex_z, slicex,'FaceColor','texturemap', 'FaceAlpha','texturemap','EdgeColor','none','AlphaDataMapping','direct','CDataMapping','direct','AlphaData',maskx);
    end

    for posy = 1:sizes(2)
        slicey_x = [0 (sizes(1)-1);0 (sizes(1)-1)]; slicey_y = [posy posy;posy posy]; slicey_z = [0 0;(sizes(3)-1) (sizes(3)-1)];
        slicey = squeeze(histo.Voxelvolume(:,posy,:)); slicey = im2uint8(slicey)';
        masky = slicey;
        surface(slicey_x,slicey_y,slicey_z, slicey,'FaceColor','texturemap', 'FaceAlpha','texturemap','EdgeColor','none','AlphaDataMapping','direct','CDataMapping','direct','AlphaData',masky);
    end

    for posz = 1:sizes(3)
        slicez_x = [0 (sizes(1)-1);0 (sizes(1)-1)]; slicez_y = [0 0;(sizes(2)-1) (sizes(2)-1)]; slicez_z = [posz posz;posz posz];
        slicez = squeeze(histo.Voxelvolume(:,:,posz)); slicez = im2uint8(slicez)';
        maskz = slicez;
        surface(slicez_x,slicez_y,slicez_z, slicez,'FaceColor','texturemap', 'FaceAlpha','texturemap','EdgeColor','none','AlphaDataMapping','direct','CDataMapping','direct','AlphaData',maskz);
    end

function lineButtonDownFcn(hObject, eventdata, handles)
histo=getMyData();
        % Get location mouse
        p = get(0, 'PointerLocation');
        pf = get(gcf, 'pos');
        p(1:2) = p(1:2)-pf(1:2);
        set(gcf, 'CurrentPoint', p(1:2));
        p = get(gca, 'CurrentPoint');
        
        % New point on mouse location
        newvalue=p(1,1)/histo.maxx; 
        
        % List for the new markers
        newvalues=zeros(1,length(histo.value)+1);
        newalphas=zeros(1,length(histo.alpha)+1);
        newcolors=zeros(length(histo.colors)+1,3);

        % Check if the new point is between old points
        index_down=find(histo.value<=newvalue); 
        if(isempty(index_down)) 
        else
            index_down=index_down(end);
            index_up=find(histo.value>newvalue); 
            if(isempty(index_up)) 
            else
                index_up=index_up(1);
                
                % Copy the (first) old markers to the new lists
                newvalues(1:index_down)=histo.value(1:index_down);
                newalphas(1:index_down)=histo.alpha(1:index_down);
                newcolors(1:index_down,:)=histo.colors(1:index_down,:);
                
                % Add the new interpolated marker
                perc=(newvalue-histo.value(index_down)) / (histo.value(index_up) - histo.value(index_down));
                color=(1-perc)*histo.colors(index_down,:)+perc*histo.colors(index_up,:);
                alpha=(1-perc)*histo.alpha(index_down)+perc*histo.alpha(index_up);
                newvalues(index_up)=newvalue; newalphas(index_up)=alpha; newcolors(index_up,:)=color;
              
                % Copy the (last) old markers to the new lists
                newvalues(index_up+1:end)=histo.value(index_up:end);
                newalphas(index_up+1:end)=histo.alpha(index_up:end);
                newcolors(index_up+1:end,:)=histo.colors(index_up:end,:);
        
                % Make the new lists the used marker lists
                histo.value=newvalues; histo.alpha=newalphas; histo.colors=newcolors;
            end
        end
        
        % Update the histogram window
        cla(histo.histogramaxeshandle);
        histo=createHistogram(histo);
        histo=drawPoints(histo,handles);
setMyData(histo);       
        
% --- Executes on mousedown on max treshold line
function pointButtonDownFcn(hObject, eventdata, handles)
histo=getMyData();
    p_sel=find(histo.pointhandle==gcbo);
    set(gcbo, 'MarkerSize',8);
    userdata=get(gcf, 'userdata');
    while userdata{1};
        % Get location mouse
        p = get(0, 'PointerLocation');
        pf = get(gcf, 'pos');
        p(1:2) = p(1:2)-pf(1:2);
        set(gcf, 'CurrentPoint', p(1:2));
        p = get(gca, 'CurrentPoint');
        
        % Set point to location mouse
        histo.value(p_sel)=p(1,1)/histo.maxx; 
        histo.alpha(p_sel)=p(1,2)/histo.maxy;
        
        % Correct new location
        if(histo.alpha(p_sel)<0), histo.alpha(p_sel)=0; end
        if(histo.alpha(p_sel)>1), histo.alpha(p_sel)=1; end
        if(histo.value(p_sel)<0), histo.value(p_sel)=0; end
        if(histo.value(p_sel)>1), histo.value(p_sel)=1; end
        if((p_sel>1)&&(histo.value(p_sel-1)>histo.value(p_sel)))
            histo.value(p_sel)=histo.value(p_sel-1);
        end
        
        if((p_sel<length(histo.value))&&(histo.value(p_sel+1)<histo.value(p_sel)))
            histo.value(p_sel)=histo.value(p_sel+1);
        end

        p = [histo.value(p_sel)*histo.maxx histo.alpha(p_sel)*histo.maxy];
        % Move point
        set(gcbo, 'xdata', p(1, 1));
        set(gcbo, 'ydata', p(1, 2));
        
        % Move line
        set(histo.linehandle, 'xdata',histo.value*histo.maxx);
        set(histo.linehandle, 'ydata',histo.alpha*histo.maxy);
        pause(.01)
        userdata=get(gcf, 'userdata');
    end
    set(hObject, 'MarkerSize', 6);

    % Create the color and alpha map
    histo=createAlphaColorMap(histo);
    
    % Update the color and alpha map
    colormap(histo.mainwindowaxeshandle,histo.ColorMap); 
    alphamap(histo.mainwindowhandle,histo.AlphaMap);
setMyData(histo);

function histo=drawPoints(histo,handles)
    % Delete old points and line
    try delete(histo.linehandle), for i=1:length(histo.pointhandle), delete(histo.pointhandle(i)), end, catch end
    
    % Display the markers and line through the markers.
    axes(histo.histogramaxeshandle);
    histo.linehandle=plot(histo.value*histo.maxx,histo.alpha*histo.maxy,'m');
    set(histo.linehandle,'ButtonDownFcn','VolumeRender(''lineButtonDownFcn'',gcbo,[],guidata(gcbo))');
    for i=1:length(histo.value)
        histo.pointhandle(i)=plot(histo.value(i)*histo.maxx,histo.alpha(i)*histo.maxy,'bo','MarkerFaceColor',histo.colors(i,:));
        set(histo.pointhandle(i),'ButtonDownFcn','VolumeRender(''pointButtonDownFcn'',gcbo,[],guidata(gcbo))');
    end
    
function data=getMyData()
    % Get the data from the GUI
    data=getappdata(gcf,'histogramdata');
    
function setMyData(data)
    % Save the data to the GUI
    setappdata(data.histogramwindowhandle,'histogramdata',data);


% --- Executes on mouse motion over figure - except title and menu.
function figure1_WindowButtonMotionFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on selection change in popupmenu_colors.
function popupmenu_colors_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_colors (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu_colors contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_colors
histo=getMyData();
    % Generate the new color markers
    c_choice=get(handles.popupmenu_colors,'Value');
    ncolors=length(histo.value);
    switch c_choice,
        case 1,new_colormap=jet(200); 
        case 2, new_colormap=hsv(200);
        case 3, new_colormap=hot(200);
        case 4, new_colormap=cool(200);
        case 5, new_colormap=spring(200);
        case 6, new_colormap=summer(200);
        case 7, new_colormap=autumn(200);
        case 8, new_colormap=winter(200);
        case 9, new_colormap=gray(200);
        case 10, new_colormap=bone(200);
        case 11, new_colormap=copper(200);
        case 12, new_colormap=pink(200);
        otherwise, new_colormap=hot(200);
    end
    new_colormap=new_colormap(round(1:(end-1)/(ncolors-1):end),:);
    histo.colors=new_colormap;
    
    % Draw the new color markers and make the color and alpha map
    histo=drawPoints(histo,handles);
    histo=createAlphaColorMap(histo);
    
    % Apply the new color and alpha map to the 3d rendering
    colormap(histo.mainwindowaxeshandle,histo.ColorMap); 
    alphamap(histo.mainwindowhandle,histo.AlphaMap);
setMyData(histo);

% --- Executes during object creation, after setting all properties.
function popupmenu_colors_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_colors (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
    set(hObject,'String',{'jet','hsv','hot','cool','spring','summer','autumn','winter','gray','bone','copper','pink'});
    set(hObject,'Value',3);
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end




