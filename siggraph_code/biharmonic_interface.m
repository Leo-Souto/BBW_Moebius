function varargout = biharmonic_interface(varargin)
% BIHARMONIC_INTERFACE MATLAB code for biharmonic_interface.fig
%      BIHARMONIC_INTERFACE, by itself, creates a new BIHARMONIC_INTERFACE or raises the existing
%      singleton*.
%
%      H = BIHARMONIC_INTERFACE returns the handle to a new BIHARMONIC_INTERFACE or the handle to
%      the existing singleton*.
%
%      BIHARMONIC_INTERFACE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BIHARMONIC_INTERFACE.M with the given input arguments.
%
%      BIHARMONIC_INTERFACE('Property','Value',...) creates a new BIHARMONIC_INTERFACE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before biharmonic_interface_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to biharmonic_interface_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help biharmonic_interface

% Last Modified by GUIDE v2.5 20-Aug-2018 23:33:09

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @biharmonic_interface_OpeningFcn, ...
    'gui_OutputFcn',  @biharmonic_interface_OutputFcn, ...
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


% --- Executes just before biharmonic_interface is made visible.
function biharmonic_interface_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to biharmonic_interface (see VARARGIN)

% Choose default command line output for biharmonic_interface

% open the image
handles.output = hObject;
handles.inputfilename = 'lamppost.jpg';
handles.matlabImage = imread(handles.inputfilename);
handles.imagesize = size(handles.matlabImage);
axes(handles.image)
handles.img_handle = image(handles.matlabImage);
axis off
% initiate the point
handles.point = [0 0];
handles.point2 = [0 0];
%collect the handles of previous window
%parece uma boa computar previamente a malha esférica
TR = IcosahedronMesh;
TR = SubdivideSphericalMesh(TR, 3);
handles.V = TR.Points;
handles.V_equi = sphere2equi(TR.Points);
handles.F = convhull(handles.V);
handles.F_equi = delaunay(handles.V_equi);

handles.num_of_handles = 0;
handles.all_handles = {};
handles.handles_plot = {};
handles.new_handle_order = {};

handles.button_value_handle = 'point_handle';
handles.bone_flag = 1;
handles.cage_flag = 0;
handles.cage_tag_counter = 0;
handles.plot_cage_tag_counter = 0;
handles.all_cages = {};
handles.all_regions = {};
handles.cross_end = [];
handles.region_positions = [];
handles.rect_flag = 0;
set(handles.figure1, 'Position', get(0, 'Screensize'));
set(handles.img_handle,'ButtonDownFcn',{@image_ButtonDownFcn,handles});
set(handles.image,'ButtonDownFcn',{@image_ButtonDownFcn,handles});

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes biharmonic_interface wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = biharmonic_interface_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in add_point_handle.
function handles = add_point_handle_Callback(hObject, eventdata, handles)
% hObject    handle to add_point_handle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

axes(handles.image)
hold on
h = plot(handles.point(1),handles.point(2),'o','MarkerFaceColor',[0.725 0.725 0.725],'MarkerEdgeColor',[0.25 0.25 0.25],'MarkerSize',7,'Parent',handles.image);
hold off
handles.num_of_handles = handles.num_of_handles + 1;
handles.handles_plot{handles.num_of_handles} = h;
equi_coord = pixel2equi(handles.point,handles.imagesize);
e_coord = sphere2equi(handles.V(knnsearch(handles.V,equi2sphere(equi_coord)),:));
%não gosto dessa parte aqui
equi_coord(2,:) = [equi_coord(1,1)+0.3 equi_coord(1,2)];
equi_coord(3,:) = [equi_coord(1,1)-0.3 equi_coord(1,2)];

handles.all_handles{end+1} = bm_handle(handles.num_of_handles,equi_coord,e_coord,'Point',0);
set(get(handles.image,'Children'),'HitTest','on','ButtonDownFcn',{@image_ButtonDownFcn,handles})
set(handles.img_handle,'ButtonDownFcn',{@image_ButtonDownFcn,handles});
% Update handles structure
guidata(hObject, handles);



% --- Executes on button press in add_curved_bone.
function handles = add_curved_bone_Callback(hObject, eventdata, handles)
% hObject    handle to add_curved_bone (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%
axes(handles.image)

equi_coord = pixel2equi(handles.point,handles.imagesize);
equi_coord2 = pixel2equi(handles.point2,handles.imagesize);
sphere_coord = equi2sphere(equi_coord);
sphere_coord2 = equi2sphere(equi_coord2);
e_coord = sphere2equi(handles.V(knnsearch(handles.V,equi2sphere(equi_coord)),:));
e_coord2 = sphere2equi(handles.V(knnsearch(handles.V,equi2sphere(equi_coord2)),:));
t = linspace(0,1,100);
parametric = [];

parametric(:,1) = (1-t)*sphere_coord(1) + t*sphere_coord2(1);
parametric(:,2) = (1-t)*sphere_coord(2) + t*sphere_coord2(2);
parametric(:,3) = (1-t)*sphere_coord(3) + t*sphere_coord2(3);
norma = sqrt(parametric(:,1).^2 + parametric(:,2).^2 + parametric(:,3).^2);
parametric(:,1) = parametric(:,1)./norma;
parametric(:,2) = parametric(:,2)./norma;
parametric(:,3) = parametric(:,3)./norma;

equi_coord_par = sphere2equi(parametric);
deriv = diff(equi_coord_par(:,1));
index = find(abs(deriv) > 0.1);
index = sort(index);

handles.num_of_handles = handles.num_of_handles + 1;
equi_coord_par_pix = equi2pixels(sphere2equi(parametric),handles.imagesize);

axes(handles.image)
hold on
if ~isempty(index)
    inicio = index(1);
    fim = index(length(index))+1;
    X1 = equi_coord_par_pix(1:inicio,1);
    Y1 = equi_coord_par_pix(1:inicio,2);
    X2 = equi_coord_par_pix(fim:100,1);
    Y2 = equi_coord_par_pix(fim:100,2);
    f = plot(X1,Y1,'Color',[0.79 0 0.125],'Parent', handles.image);
    f2 = plot(X2,Y2,'Color',[0.79 0 0.125],'Parent', handles.image);
    handles.handles_plot{end+1} = [f f2];
    handles.cross_end = [handles.cross_end handles.num_of_handles];
else
    f = plot(equi_coord_par_pix(:,1),equi_coord_par_pix(:,2),'Color',[0.79 0 0.125]);
    handles.handles_plot{end+1} = f;
end
h =  plot(handles.point2(1),handles.point2(2),'o','MarkerFaceColor',[0.79 0 0.125],'MarkerEdgeColor',[0.25 0.25 0.25],'MarkerSize',5,'Parent',handles.image);
hold off
handles.handles_plot{end+1} = h;
if isempty(index)
    for k = 1:100
        prop = bm_length_of_curve(equi_coord_par,k)/bm_length_of_curve(equi_coord_par,100);
        if prop >= 0.5
            ind = k;
            break
        end
    end
    handles.all_handles{handles.num_of_handles} = bm_handle(handles.num_of_handles,[equi_coord;...
        equi_coord_par(ind,1) equi_coord_par(ind,2);equi_coord2],[e_coord;e_coord2],'Curved',0);
    handles.all_handles{handles.num_of_handles}.proportion = bm_length_of_curve(equi_coord_par,ind)/bm_length_of_curve(equi_coord_par,100);
else
    ind = 0;
    A = length(equi_coord_par(1:inicio,:));
    B = length(equi_coord_par(fim:end,:));
    for k = 1:inicio
        prop = bm_length_of_curve(equi_coord_par,k)/...
            (bm_length_of_curve(equi_coord_par,A)+ bm_length_of_curve(equi_coord_par(fim:end,:),B));
        if prop >= 0.5
            ind = k;
            break
        end
    end
    if ind == 0
        ind = floor(4*inicio/5);
    end
    handles.all_handles{handles.num_of_handles} = bm_handle(handles.num_of_handles,[equi_coord;...
        equi_coord_par(ind,1) equi_coord_par(ind,2);equi_coord2],[e_coord;e_coord2],'Curved',0);
    handles.all_handles{handles.num_of_handles}.proportion = bm_length_of_curve(equi_coord_par,ind)/...
        (bm_length_of_curve(equi_coord_par,A)+ bm_length_of_curve(equi_coord_par(fim:end,:),B));
end
set(get(handles.image,'Children'),'HitTest','on','ButtonDownFcn',{@image_ButtonDownFcn,handles})
set(handles.img_handle,'ButtonDownFcn',{@image_ButtonDownFcn,handles},'HitTest','on');
guidata(hObject,handles)

function handles = add_cage_handle_Callback(hObject,eventdata,handles)

if handles.cage_flag == 0
    axes(handles.image)
    hold on
    h = plot(handles.point(1),handles.point(2),'bo','MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[0.25 0.25 0.25],'MarkerSize',5,'Parent',handles.image);
    hold off
    handles.num_of_handles = handles.num_of_handles + 1;
    handles.handles_plot{handles.num_of_handles} = h;
    equi_coord = pixel2equi(handles.point,handles.imagesize);
    equi_coord = [reshape(equi_coord,[1 2]); equi_coord(1)+0.3 equi_coord(2); equi_coord(1)-0.3 equi_coord(2)];
    e_coord = sphere2equi(handles.V(knnsearch(handles.V,equi2sphere(equi_coord)),:));
    if isequal(handles.button_value_handle,'closed_cage_handle')
        handles.all_handles{end+1} = bm_handle(handles.num_of_handles,equi_coord,e_coord,'Closed_Cage',1);
    end
    handles.cage_tag_counter = handles.cage_tag_counter + 1;
    handles.all_cages{end+1} = handles.num_of_handles;
    handles.plot_cage_tag_counter = 1;
    handles.all_handles{end}.cage_tag = handles.cage_tag_counter;
    set(get(handles.image,'Children'),'HitTest','on','ButtonDownFcn',{@image_ButtonDownFcn,handles})
    set(handles.img_handle,'ButtonDownFcn',{@image_ButtonDownFcn,handles});
else
    
    
    equi_coord = pixel2equi(handles.point,handles.imagesize);
    equi_coord2 = pixel2equi(handles.point2,handles.imagesize);
    sphere_coord = equi2sphere(equi_coord);
    sphere_coord2 = equi2sphere(equi_coord2);
    t = linspace(0,1,100);
    parametric = [];
    
    parametric(:,1) = (1-t)*sphere_coord(1) + t*sphere_coord2(1);
    parametric(:,2) = (1-t)*sphere_coord(2) + t*sphere_coord2(2);
    parametric(:,3) = (1-t)*sphere_coord(3) + t*sphere_coord2(3);
    norma = sqrt(parametric(:,1).^2 + parametric(:,2).^2 + parametric(:,3).^2);
    parametric(:,1) = parametric(:,1)./norma;
    parametric(:,2) = parametric(:,2)./norma;
    parametric(:,3) = parametric(:,3)./norma;
    
    equi_coord_par = sphere2equi(parametric);
    deriv = diff(equi_coord_par(:,1));
    index = find(abs(deriv) > 0.1);
    index = sort(index);
    equi_coord_par_pix = equi2pixels(sphere2equi(parametric),handles.imagesize);
    axes(handles.image)
    hold on
    if ~isempty(index)
        inicio = index(1);
        fim = index(length(index))+1;
        X1 = equi_coord_par_pix(1:inicio,1);
        Y1 = equi_coord_par_pix(1:inicio,2);
        X2 = equi_coord_par_pix(fim:100,1);
        Y2 = equi_coord_par_pix(fim:100,2);
        f = plot(X1,Y1,'Color',[1 1 1],'Parent', handles.image);
        f2 = plot(X2,Y2,'Color',[1 1 1],'Parent', handles.image);
        handles.handles_plot{end+1} = [f f2];
        handles.cross_end = [handles.cross_end handles.num_of_handles];
    else
        f = line('Parent', handles.image, 'XData', equi_coord_par_pix(:,1),...
            'YData', equi_coord_par_pix(:,2),'Color',[1 1 1]);
        handles.handles_plot{end+1} = f;
    end
    h = plot(handles.point2(1),handles.point2(2),'bo','MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[0.25 0.25 0.25],'MarkerSize',5,'Parent',handles.image);
    handles.handles_plot{end} = [handles.handles_plot{end} h];
    hold off
    handles.num_of_handles = handles.num_of_handles + 1;
    equi_coord = pixel2equi(handles.point2,handles.imagesize);
    equi_coord = [reshape(equi_coord,[1 2]); equi_coord(1)+0.3 equi_coord(2); equi_coord(1)-0.3 equi_coord(2)];
    e_coord = sphere2equi(handles.V(knnsearch(handles.V,equi2sphere(equi_coord)),:));
    handles.plot_cage_tag_counter = handles.plot_cage_tag_counter + 1;
    if isequal(handles.button_value_handle,'closed_cage_handle')
        handles.all_handles{end+1} = bm_handle(handles.num_of_handles,equi_coord,e_coord,'Closed_Cage',handles.plot_cage_tag_counter);
    end
    handles.all_handles{end}.cage_tag = handles.cage_tag_counter;
    handles.all_cages{handles.cage_tag_counter} = [handles.all_cages{handles.cage_tag_counter} handles.num_of_handles];
    
end


% --- Executes on mouse press over axes background.
function image_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to image (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pt = get(get(hObject,'Parent'),'CurrentPoint');

if isequal(handles.button_value_handle,'point_handle') && isequal(get(gcf,'SelectionType'),'normal')
    handles.point = [pt(1,1) pt(1,2)];
    handles = add_point_handle_Callback(hObject,eventdata,handles);
elseif isequal(handles.button_value_handle,'bone_handle') && handles.bone_flag == 1 && isequal(get(gcf,'SelectionType'),'normal')
    handles.bone_flag = 2;
    handles.point = [pt(1,1) pt(1,2)];
    axes(handles.image)
    hold on
    g =  plot(handles.point(1),handles.point(2),'o','MarkerFaceColor',[0.79 0 0.125],'MarkerEdgeColor',[0.25 0.25 0.25],'MarkerSize',5,'Parent',handles.image);
    handles.handles_plot{end+1} = g;
    hold off
elseif isequal(handles.button_value_handle,'bone_handle') && handles.bone_flag == 2 && isequal(get(gcf,'SelectionType'),'normal')
    handles.bone_flag = 1;
    handles.point2 = [pt(1,1) pt(1,2)];
    handles = add_curved_bone_Callback(hObject,eventdata,handles);
elseif isequal(handles.button_value_handle,'closed_cage_handle') && isequal(get(gcf,'SelectionType'),'normal')
    if handles.cage_flag == 0
        handles.point = [pt(1,1) pt(1,2)];
        handles.point2 = handles.point;
    else
        handles.point = handles.point2;
        handles.point2 = [pt(1,1) pt(1,2)];
    end
    handles = add_cage_handle_Callback(hObject,eventdata,handles);
    handles.cage_flag = 1;
elseif isequal(handles.button_value_handle,'conformal_region')
    handles.rect{1} = imrect;
    pos = getPosition(handles.rect{1});
    x_pos = [pos(1) pos(1) pos(1)+pos(3) pos(1)+pos(3)];
    y_pos = [pos(2) pos(2)+pos(4) pos(2)+pos(4) pos(2)];
    handles.all_regions{end+1} = patch(x_pos,y_pos,'g',...
        'FaceColor', [0.25 0.25 0.25], 'FaceAlpha', 0.4,'EdgeColor',[0.725 0.725 0.725],'Parent',handles.image);
    delete(handles.rect{1})
    handles.region_positions = [handles.region_positions; pos];
    
end
set(get(handles.image,'Children'),'HitTest','on','ButtonDownFcn',{@image_ButtonDownFcn,handles})
set(handles.img_handle,'ButtonDownFcn',{@image_ButtonDownFcn,handles});
guidata(hObject,handles);

% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function figure1_WindowButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isequal(handles.button_value_handle,'closed_cage_handle') && isequal(get(gcf,'SelectionType'),'alt') && ...
        length(handles.all_cages{handles.cage_tag_counter}) > 2
    index = handles.all_cages{length(handles.all_cages)}(1);
    first_handle = handles.all_handles{index};
    handles.point = equi2pixels(first_handle.position(1,:),handles.imagesize);
    equi_coord = pixel2equi(handles.point,handles.imagesize);
    equi_coord2 = pixel2equi(handles.point2,handles.imagesize);
    sphere_coord = equi2sphere(equi_coord);
    sphere_coord2 = equi2sphere(equi_coord2);
    t = linspace(0,1,100);
    parametric = [];
    
    parametric(:,1) = (1-t)*sphere_coord(1) + t*sphere_coord2(1);
    parametric(:,2) = (1-t)*sphere_coord(2) + t*sphere_coord2(2);
    parametric(:,3) = (1-t)*sphere_coord(3) + t*sphere_coord2(3);
    norma = sqrt(parametric(:,1).^2 + parametric(:,2).^2 + parametric(:,3).^2);
    parametric(:,1) = parametric(:,1)./norma;
    parametric(:,2) = parametric(:,2)./norma;
    parametric(:,3) = parametric(:,3)./norma;
    
    equi_coord_par = sphere2equi(parametric);
    deriv = diff(equi_coord_par(:,1));
    index = find(abs(deriv) > 0.1);
    index = sort(index);
    equi_coord_par_pix = equi2pixels(sphere2equi(parametric),handles.imagesize);
    axes(handles.image)
    hold on
    if ~isempty(index)
        inicio = index(1);
        fim = index(length(index))+1;
        X1 = equi_coord_par_pix(1:inicio,1);
        Y1 = equi_coord_par_pix(1:inicio,2);
        X2 = equi_coord_par_pix(fim:100,1);
        Y2 = equi_coord_par_pix(fim:100,2);
        f = plot(X1,Y1,'Color',[1 1 1],'Parent', handles.image);
        f2 = plot(X2,Y2,'Color',[1 1 1],'Parent', handles.image);
        handles.handles_plot{end+1} = [f f2];
        %         handles.cross_end = [handles.cross_end handles.num_of_handles];
    else
        f = line('Parent', handles.image, 'XData', equi_coord_par_pix(:,1),...
            'YData', equi_coord_par_pix(:,2),'Color',[1 1 1]);
            handles.handles_plot{end+1} = f;
    end
    hold off
    handles.all_cages{handles.cage_tag_counter} = [handles.all_cages{handles.cage_tag_counter}, handles.all_cages{length(handles.all_cages)}(1)];
    handles.closed_cage_handle = 0;
    handles.cage_flag = 0;
elseif isequal(handles.button_value_handle,'closed_cage_handle') && isequal(get(gcf,'SelectionType'),'alt') && ...
        length(handles.all_cages{handles.cage_tag_counter}) == 2
    handles.closed_cage_handle = 0;
    handles.cage_flag = 0;
end

set(handles.img_handle,'ButtonDownFcn',{@image_ButtonDownFcn,handles});
set(handles.image,'ButtonDownFcn',{@image_ButtonDownFcn,handles});
guidata(hObject,handles);

% --- Executes when selected object is changed in uibuttongroup1.
function uibuttongroup1_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uibuttongroup1
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%if button value handle == closed cage atualizar o fim da cage

handles.button_value_handle = get(eventdata.NewValue, 'Tag');
handles.closed_cage_handle = 0;
handles.cage_flag = 0;
handles.bone_flag = 1;

set(get(handles.image,'Children'),'HitTest','on','ButtonDownFcn',{@image_ButtonDownFcn,handles})
set(handles.img_handle,'ButtonDownFcn',{@image_ButtonDownFcn,handles});
set(handles.image,'ButtonDownFcn',{@image_ButtonDownFcn,handles});
guidata(hObject,handles);



% --- Executes on button press in delete_handle.
function delete_handle_Callback(hObject, eventdata, handles)
% hObject    handle to delete_handle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.image)
children = get(gca, 'children');
if handles.num_of_handles > 0
    switch handles.all_handles{handles.num_of_handles}.type
        case 'Curved'
            if isempty(find(handles.cross_end == handles.num_of_handles)) && handles.bone_flag == 1
                delete(children(1:3))
                handles.handles_plot(end) = [];
                handles.handles_plot(end) = [];
                handles.all_handles(end) = [];
                handles.num_of_handles = handles.num_of_handles - 1;
            elseif ~isempty(find(handles.cross_end == handles.num_of_handles)) && handles.bone_flag == 1
                delete(children(1:4))
                handles.cross_end = handles.cross_end(handles.cross_end ~= handles.num_of_handles);
                handles.handles_plot(end) = [];
                handles.handles_plot(end) = [];
                handles.all_handles(end) = [];
                handles.num_of_handles = handles.num_of_handles - 1;
            elseif handles.bone_flag == 2
                handles.bone_flag = 1;
                delete(children(1))
            end
        case 'Point'
            delete(children(1))
            handles.handles_plot(end) = [];
            handles.all_handles(end) = [];
            handles.num_of_handles = handles.num_of_handles - 1;
        case 'Closed_Cage'
            if handles.cage_flag == 0
                delete(children(1:3))
                handles.cage_flag = 1;
                handles.all_cages{handles.cage_tag_counter} = handles.all_cages{handles.cage_tag_counter}(1:end-2);
                handles.closed_cage_handle = length(handles.all_cages{handles.cage_tag_counter});
                handles.plot_cage_tag_counter = handles.plot_cage_tag_counter - 1;
                handles.handles_plot(end) = [];
                handles.all_handles(end) = [];
                handles.num_of_handles = handles.num_of_handles - 1;
            else
                if length(handles.all_cages{handles.cage_tag_counter}) > 1
                    delete(children(1:2))
                    handles.all_cages{handles.cage_tag_counter} = handles.all_cages{handles.cage_tag_counter}(1:end-1);
                    handles.closed_cage_handle = length(handles.all_cages{handles.cage_tag_counter});
                    handles.plot_cage_tag_counter = handles.plot_cage_tag_counter - 1;
                    handles.handles_plot(end) = [];
                    handles.all_handles(end) = [];
                    handles.num_of_handles = handles.num_of_handles - 1;
                else
                    delete(children(1))
                    handles.cage_tag_counter = handles.cage_tag_counter - 1;
                    handles.cage_flag = 0;
                    handles.all_cages(end) = [];
                    handles.handles_plot(end) = [];
                    handles.all_handles(end) = [];
                    handles.num_of_handles = handles.num_of_handles - 1;
                end
                
            end
            if handles.num_of_handles > 2
                handles.point2 = equi2pixels(handles.all_handles{handles.num_of_handles - 1}.position(1,:),handles.imagesize);
                handles.point = equi2pixels(handles.all_handles{handles.num_of_handles - 2}.position(1,:),handles.imagesize);
                handles.handles_plot(end) = [];
                handles.all_handles(end) = [];
                handles.num_of_handles = handles.num_of_handles - 1;
            elseif handles.num_of_handles == 2
                handles.point2 = equi2pixels(handles.all_handles{handles.num_of_handles - 1}.position(1,:),handles.imagesize);
                handles.point = handles.point2;
                handles.handles_plot(end) = [];
                handles.all_handles(end) = [];
                handles.num_of_handles = handles.num_of_handles - 1;
            end
    end
elseif handles.num_of_handles == 0 && isequal(handles.button_value_handle,'bone_handle') && handles.bone_flag == 2
    handles.bone_flag = 1;
    delete(children(1))
end

setappdata(0,'Handles',handles);
guidata(hObject,handles)


% --- Executes on button press in next.
function next_Callback(hObject, eventdata, handles)
% hObject    handle to next (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
positions = [];
% h = warndlg('The next window may take a while to open due to weights optmization.','Warning');
% uiwait(h);

%% Here we put the handles in order [Points, then bones]
num = handles.num_of_handles;
for i = 1:num
    positions = cat(1,positions,handles.all_handles{i}.mesh_position);
end
point_indices = [];
curved_bone_indices = [];
cage_indices = [];
new_indice_order_cages = [];
for i = 1:num
    switch handles.all_handles{i}.type
        case 'Point'
            indice = knnsearch(positions,handles.all_handles{i}.mesh_position(1,:));
            point_indices = [point_indices,indice];
            handles.new_handle_order{end+1} = handles.all_handles{i};
        case 'Closed_Cage'
            indice = knnsearch(positions,handles.all_handles{i}.mesh_position(1,:));
            point_indices = [point_indices,indice];
            cage_indices = [cage_indices; i length(point_indices)];
            handles.new_handle_order{end+1} = handles.all_handles{i};
            new_indice_order_cages = [new_indice_order_cages; i length(handles.new_handle_order)];
    end
end
for i = 1:num
    switch handles.all_handles{i}.type
        case 'Curved'
            indice1 = knnsearch(positions,handles.all_handles{i}.mesh_position(1,:));
            indice2 = knnsearch(positions,handles.all_handles{i}.mesh_position(2,:));
            curved_bone_indices = cat(1,curved_bone_indices,[indice1 indice2]);
            handles.new_handle_order{end+1} = handles.all_handles{i};
    end
end
not_finally_cage_indices = cell(1, length(handles.all_cages));
finally_cage_indices = [];
handles.new_cages = cell(1, length(handles.all_cages));
if ~isempty(handles.all_cages)
    for i = 1:length(handles.all_cages)
        for j = handles.all_cages{i}
            not_finally_cage_indices{i} = [not_finally_cage_indices{i} cage_indices((find(cage_indices(:,1) == j)),2)];
            handles.new_cages{i} = [handles.new_cages{i} new_indice_order_cages(find(new_indice_order_cages(:,1) == j),2)];
        end
    end
    for i = 1:length(handles.all_cages)
        for j = 1:length(not_finally_cage_indices{i})-1
            finally_cage_indices = [finally_cage_indices; not_finally_cage_indices{i}(j) not_finally_cage_indices{i}(j+1)];
        end
    end
end

%% Here we take care of nearby endpoints
handles.coupled_endpoints = cell(handles.num_of_handles,1);
handles.handles_visibility = ones(handles.num_of_handles,3);
handles.bone_intersections = zeros(handles.num_of_handles,handles.num_of_handles);
handles.intersection_points = [];
handles.point_coord_clusters = [];
handles.point_index_clusters = [];

for i  = 1:handles.num_of_handles
    if isequal(handles.new_handle_order{i}.type,'Curved')
        P1 = handles.new_handle_order{i}.position(1,:);
        P3 = handles.new_handle_order{i}.position(3,:);
        handles.point_coord_clusters = [handles.point_coord_clusters; P1 ; P3];
        handles.point_index_clusters = [handles.point_index_clusters; i 1; i 3];
    end
end
raio_approx = 2*pi/100;
[clustersCentroids,clustersGeoMedians,clustersXY, handles.clustersInd] = clusterXYpoints(handles.point_coord_clusters,raio_approx,1,'centroid','merge');
clustersCentroidsMeshPos = sphere2equi(handles.V(knnsearch(handles.V,equi2sphere(clustersCentroids)),:));
handles.num_of_clusters = numel(handles.clustersInd);
for i = 1:handles.num_of_clusters
    for j = handles.clustersInd{i}
        index = handles.point_index_clusters(j,:);
        handles.new_handle_order{index(1)}.position(index(2),:) = clustersCentroids(i,:);
        handles.new_handle_order{index(1)}.new_position = handles.new_handle_order{index(1)}.position;
        handles.new_handle_order{index(1)}.mesh_position(index(2),:) = clustersCentroidsMeshPos(i,:);
    end
end
for i = 1:handles.num_of_clusters
    if length(handles.clustersInd{i}) > 1
        for j = handles.clustersInd{i}(2:end)
            index = handles.point_index_clusters(j,:);
            handles.handles_visibility(index(1),index(2)) = 0;
        end
    end
end
for i = 1:length(handles.new_handle_order)
    if isequal(handles.new_handle_order{i}.type,'Curved')
        if isequal(handles.new_handle_order{i}.mesh_position(1,:),handles.new_handle_order{i}.mesh_position(3,:))
            handles.new_handle_order{i} = {};
        end
    end
end
handles.new_handle_order = handles.new_handle_order(~cellfun('isempty',handles.new_handle_order));
handles.num_of_handles = length(handles.new_handle_order);
%% Here we take care of nearby endpoints
handles.coupled_endpoints = cell(handles.num_of_handles,1);
handles.handles_visibility = ones(handles.num_of_handles,3);
handles.bone_intersections = zeros(handles.num_of_handles,handles.num_of_handles);
handles.intersection_points = [];
handles.point_coord_clusters = [];
handles.point_index_clusters = [];

for i  = 1:handles.num_of_handles
    if isequal(handles.new_handle_order{i}.type,'Curved')
        P1 = handles.new_handle_order{i}.position(1,:);
        P3 = handles.new_handle_order{i}.position(3,:);
        handles.point_coord_clusters = [handles.point_coord_clusters; P1 ; P3];
        handles.point_index_clusters = [handles.point_index_clusters; i 1; i 3];
    end
end
raio_approx = 2*pi/100;
[clustersCentroids,clustersGeoMedians,clustersXY, handles.clustersInd] = clusterXYpoints(handles.point_coord_clusters,raio_approx,1,'centroid','merge');
clustersCentroidsMeshPos = sphere2equi(handles.V(knnsearch(handles.V,equi2sphere(clustersCentroids)),:));
handles.num_of_clusters = numel(handles.clustersInd);
for i = 1:handles.num_of_clusters
    for j = handles.clustersInd{i}
        index = handles.point_index_clusters(j,:);
        handles.new_handle_order{index(1)}.position(index(2),:) = clustersCentroids(i,:);
        handles.new_handle_order{index(1)}.new_position = handles.new_handle_order{index(1)}.position;
        handles.new_handle_order{index(1)}.mesh_position(index(2),:) = clustersCentroidsMeshPos(i,:);
    end
end
for i = 1:handles.num_of_clusters
    if length(handles.clustersInd{i}) > 1
        for j = handles.clustersInd{i}(2:end)
            index = handles.point_index_clusters(j,:);
            handles.handles_visibility(index(1),index(2)) = 0;
        end
    end
end
num_of_intersections = cell(handles.num_of_handles,1);
%%[x_int y_int ind_linspace]

%%in this desnecessary big code we break the bones
for i =1:handles.num_of_handles
    if strcmp(handles.new_handle_order{i}.type,'Curved')
        for j = 1:handles.num_of_handles
            if strcmp(handles.new_handle_order{j}.type,'Curved') && i ~= j
                
                segmento_i_pos = [handles.new_handle_order{i}.position(1,:); ...
                    handles.new_handle_order{i}.position(3,:)];
                segmento_j_pos = [handles.new_handle_order{j}.position(1,:); ...
                    handles.new_handle_order{j}.position(3,:)];
                
                sphere_i_pos = equi2sphere(segmento_i_pos);
                sphere_j_pos = equi2sphere(segmento_j_pos);
                segmento_i_x = linspace(sphere_i_pos(1,1),sphere_i_pos(2,1),100);
                segmento_i_y = linspace(sphere_i_pos(1,2),sphere_i_pos(2,2),100);
                segmento_i_z = linspace(sphere_i_pos(1,3),sphere_i_pos(2,3),100);
                segmento_j_x = linspace(sphere_j_pos(1,1),sphere_j_pos(2,1),100);
                segmento_j_y = linspace(sphere_j_pos(1,2),sphere_j_pos(2,2),100);
                segmento_j_z = linspace(sphere_j_pos(1,3),sphere_j_pos(2,3),100);
                segmento_i = [segmento_i_x' segmento_i_y' segmento_i_z'];
                segmento_j = [segmento_j_x' segmento_j_y' segmento_j_z'];
                segmento_i_norm = sqrt(sum(segmento_i.^2,2));
                segmento_i_norm = [segmento_i_norm segmento_i_norm segmento_i_norm];%%faço isso pra poder usar div-by-element
                segmento_j_norm = sqrt(sum(segmento_j.^2,2));
                segmento_j_norm = [segmento_j_norm segmento_j_norm segmento_j_norm];
                segmento_i = segmento_i./segmento_i_norm;
                segmento_j = segmento_j./segmento_j_norm;
                segmento_i = sphere2equi(segmento_i);%aqui tenho o segmento reconstruido no espaço equirretangular
                segmento_j = sphere2equi(segmento_j);
                deriv_i = diff(segmento_i(:,1));
                index_i = find(abs(deriv_i) > 0.1);
                index_i = sort(index_i);
                deriv_j = diff(segmento_j(:,1));
                index_j = find(abs(deriv_j) > 0.1);
                index_j = sort(index_j);
                
                if ~isempty(index_i)
                    if ~isempty(index_j)
                        %%i and j intersect the boundary
                        seg_i_1 = segmento_i(1:index_i(1),:);
                        seg_i_2 = segmento_i(index_i(end)+1:100,:);
                        seg_j_1 = segmento_j(1:index_j(1),:);
                        seg_j_2 = segmento_j(index_j(end)+1:100,:);
                        [x1,y1] = intersections(seg_i_1(:,1),seg_i_1(:,2),seg_j_1(:,1),seg_j_1(:,2),'Robust');
                        [x2,y2] = intersections(seg_i_2(:,1),seg_i_2(:,2),seg_j_1(:,1),seg_j_1(:,2),'Robust');
                        [x3,y3] = intersections(seg_i_1(:,1),seg_i_1(:,2),seg_j_2(:,1),seg_j_2(:,2),'Robust');
                        [x4,y4] = intersections(seg_i_2(:,1),seg_i_2(:,2),seg_j_2(:,1),seg_j_2(:,2),'Robust');
                        x_intersect = [x1 x2 x3 x4];
                        y_intersect = [y1 y2 y3 y4];
                    else
                        %%only i intersect the boundary
                        seg_i_1 = segmento_i(1:index_i(1),:);
                        seg_i_2 = segmento_i(index_i(end)+1:100,:);
                        [x1,y1] = intersections(seg_i_1(:,1),seg_i_1(:,2),segmento_j(:,1),segmento_j(:,2),'Robust');
                        [x2,y2] = intersections(seg_i_2(:,1),seg_i_2(:,2),segmento_j(:,1),segmento_j(:,2),'Robust');
                        x_intersect = [x1 x2];
                        y_intersect = [y1 y2];
                    end
                else
                    if ~isempty(index_j)
                        %%only j intersect the boundary
                        seg_j_1 = segmento_j(1:index_j(1),:);
                        seg_j_2 = segmento_j(index_j(end)+1:100,:);
                        [x1,y1] = intersections(segmento_i(:,1),segmento_i(:,2),seg_j_1(:,1),seg_j_1(:,2),'Robust');
                        [x2,y2] = intersections(segmento_i(:,1),segmento_i(:,2),seg_j_2(:,1),seg_j_2(:,2),'Robust');
                        x_intersect = [x1 x2];
                        y_intersect = [y1 y2];
                    else
                        %%none intersect the boundary
                        [x_intersect,y_intersect] = intersections(segmento_i(:,1),segmento_i(:,2),segmento_j(:,1),segmento_j(:,2),'Robust');
                    end
                end
                %test if is not a simple intersection of
                %endpoints
                pos_i_1 = handles.new_handle_order{i}.position(1,:);
                pos_i_3 = handles.new_handle_order{i}.position(3,:);
                pos_j_1 = handles.new_handle_order{j}.position(1,:);
                pos_j_3 = handles.new_handle_order{j}.position(3,:);
                cond11 = 0;
                cond13 = 0;
                cond31 = 0;
                cond33 = 0;
                if ~isempty(x_intersect) && isequal(size(x_intersect),[1 1])
                    cond11 = sqrt((pos_i_1(1) - x_intersect)^2 + (pos_i_1(2) - y_intersect)^2) < 1e-5 && ...
                        sqrt((pos_j_1(1) - x_intersect)^2 + (pos_j_1(2) - y_intersect)^2) < 1e-5;
                    cond31 = sqrt((pos_i_3(1) - x_intersect)^2 + (pos_i_3(2) - y_intersect)^2) < 1e-5 && ...
                        sqrt((pos_j_1(1) - x_intersect)^2 + (pos_j_1(2) - y_intersect)^2) < 1e-5;
                    cond13 = sqrt((pos_i_1(1) - x_intersect)^2 + (pos_i_1(2) - y_intersect)^2) < 1e-5 && ...
                        sqrt((pos_j_3(1) - x_intersect)^2 + (pos_j_3(2) - y_intersect)^2) < 1e-5;
                    cond33 = sqrt((pos_i_3(1) - x_intersect)^2 + (pos_i_3(2) - y_intersect)^2) < 1e-5 && ...
                        sqrt((pos_j_3(1) - x_intersect)^2 + (pos_j_3(2) - y_intersect)^2) < 1e-5;
                end
                cond = cond11 || cond31 || cond13 || cond33;
                if ~isempty(x_intersect) && ~cond
                    ind_lin = knnsearch(segmento_i,[x_intersect y_intersect]);
                    num_of_intersections{i} = [num_of_intersections{i};...
                        x_intersect y_intersect ind_lin];
                end
            end
        end
    end
end

new_bones = 0;
for i = 1:handles.num_of_handles
    if size(num_of_intersections{i},1) > 1
        num_of_intersections{i} = sortrows(num_of_intersections{i},3);
    end
end

for i = 1:handles.num_of_handles
    if ~isempty(num_of_intersections{i})
        final_bone_pos = handles.new_handle_order{i}.position(3,:);
        init_pos = handles.new_handle_order{i}.position(1,:);
        final_pos = num_of_intersections{i}(1,1:2);
        segmento = [init_pos;final_pos];
        segmento_sph = equi2sphere(segmento);
        seg_x = linspace(segmento_sph(1,1),segmento_sph(2,1),100);
        seg_y = linspace(segmento_sph(1,2),segmento_sph(2,2),100);
        seg_z = linspace(segmento_sph(1,3),segmento_sph(2,3),100);
        seg = [seg_x' seg_y' seg_z'];
        seg_norm = sqrt(sum(seg.^2,2));
        seg_norm = [seg_norm seg_norm seg_norm];
        seg = seg./seg_norm;
        segmento = sphere2equi(seg);
        
        deriv = diff(segmento(:,1));
        index = find(abs(deriv) > 0.1);
        index = sort(index);
        if isempty(index)
            for k = 1:100
                prop = bm_length_of_curve(segmento,k)/bm_length_of_curve(segmento,100);
                if prop >= 0.5
                    ind = k;
                    break
                end
            end
            midpoint = segmento(ind,:);
            equi_coord = [init_pos;midpoint;final_pos];
            mesh_coord = sphere2equi(handles.V(knnsearch(handles.V,equi2sphere(equi_coord)),:));
            handles.new_handle_order{i} = bm_handle(handles.num_of_handles+new_bones,equi_coord,mesh_coord,'Curved',0);
            handles.new_handle_order{i}.proportion = bm_length_of_curve(segmento,ind)/bm_length_of_curve(segmento,100);
        else
            inicio = index(1);
            fim = index(length(index))+1;
            ind = 0;
            A = length(segmento(1:inicio,:));
            B = length(segmento(fim:end,:));
            for k = 1:inicio
                prop = bm_length_of_curve(segmento,k)/...
                    (bm_length_of_curve(segmento,A)+ bm_length_of_curve(segmento(fim:end,:),B));
                if prop >= 0.5
                    ind = k;
                    break
                end
            end
            if ind == 0
                ind = floor(4*inicio/5);
            end
            midpoint = segmento(ind,:);
            equi_coord = [init_pos;midpoint;final_pos];
            mesh_coord = sphere2equi(handles.V(knnsearch(handles.V,equi2sphere(equi_coord)),:));
            handles.new_handle_order{i} = bm_handle(handles.num_of_handles+new_bones,equi_coord,mesh_coord,'Curved',0);
            handles.new_handle_order{i}.proportion = bm_length_of_curve(segmento,ind)/...
                (bm_length_of_curve(segmento,A)+ bm_length_of_curve(segmento(fim:end,:),B));
        end
        if size(num_of_intersections{i},1) > 1
            for j = 2:size(num_of_intersections{i},1)
                init_pos = final_pos;
                final_pos = num_of_intersections{i}(j,1:2);
                segmento = [init_pos;final_pos];
                segmento_sph = equi2sphere(segmento);
                seg_x = linspace(segmento_sph(1,1),segmento_sph(2,1),100);
                seg_y = linspace(segmento_sph(1,2),segmento_sph(2,2),100);
                seg_z = linspace(segmento_sph(1,3),segmento_sph(2,3),100);
                seg = [seg_x' seg_y' seg_z'];
                seg_norm = sqrt(sum(seg.^2,2));
                seg_norm = [seg_norm seg_norm seg_norm];
                seg = seg./seg_norm;
                segmento = sphere2equi(seg);
                new_bones = new_bones+1;
                deriv = diff(segmento(:,1));
                index = find(abs(deriv) > 0.1);
                index = sort(index);
                if isempty(index)
                    for k = 1:100
                        prop = bm_length_of_curve(segmento,k)/bm_length_of_curve(segmento,100);
                        if prop >= 0.5
                            ind = k;
                            break
                        end
                    end
                    midpoint = segmento(ind,:);
                    equi_coord = [init_pos;midpoint;final_pos];
                    mesh_coord = sphere2equi(handles.V(knnsearch(handles.V,equi2sphere(equi_coord)),:));
                    handles.new_handle_order{end+1} = bm_handle(handles.num_of_handles+new_bones,equi_coord,mesh_coord,'Curved',0);
                    handles.new_handle_order{end}.proportion = bm_length_of_curve(segmento,ind)/bm_length_of_curve(segmento,100);
                else
                    inicio = index(1);
                    fim = index(length(index))+1;
                    ind = 0;
                    A = length(segmento(1:inicio,:));
                    B = length(segmento(fim:end,:));
                    for k = 1:inicio
                        prop = bm_length_of_curve(segmento,k)/...
                            (bm_length_of_curve(segmento,A)+ bm_length_of_curve(segmento(fim:end,:),B));
                        if prop >= 0.5
                            ind = k;
                            break
                        end
                    end
                    if ind == 0
                        ind = floor(4*inicio/5);
                    end
                    midpoint = segmento(ind,:);
                    equi_coord = [init_pos;midpoint;final_pos];
                    handles.new_handle_order{end+1} = bm_handle(handles.num_of_handles+new_bones,equi_coord,mesh_coord,'Curved',0);
                    handles.new_handle_order{end}.proportion = bm_length_of_curve(segmento,ind)/...
                        (bm_length_of_curve(segmento,A)+ bm_length_of_curve(segmento(fim:end,:),B));
                end
            end
        end
        init_pos = final_pos;
        final_pos = final_bone_pos;
        segmento = [init_pos;final_pos];
        segmento_sph = equi2sphere(segmento);
        seg_x = linspace(segmento_sph(1,1),segmento_sph(2,1),100);
        seg_y = linspace(segmento_sph(1,2),segmento_sph(2,2),100);
        seg_z = linspace(segmento_sph(1,3),segmento_sph(2,3),100);
        seg = [seg_x' seg_y' seg_z'];
        seg_norm = sqrt(sum(seg.^2,2));
        seg_norm = [seg_norm seg_norm seg_norm];
        seg = seg./seg_norm;
        segmento = sphere2equi(seg);
        new_bones = new_bones+1;
        deriv = diff(segmento(:,1));
        index = find(abs(deriv) > 0.1);
        index = sort(index);
        if isempty(index)
            for k = 1:100
                prop = bm_length_of_curve(segmento,k)/bm_length_of_curve(segmento,100);
                if prop >= 0.5
                    ind = k;
                    break
                end
            end
            midpoint = segmento(ind,:);
            equi_coord = [init_pos;midpoint;final_pos];
            mesh_coord = sphere2equi(handles.V(knnsearch(handles.V,equi2sphere(equi_coord)),:));
            handles.new_handle_order{end+1} = bm_handle(handles.num_of_handles+new_bones,equi_coord,mesh_coord,'Curved',0);
            handles.new_handle_order{end}.proportion = bm_length_of_curve(segmento,ind)/bm_length_of_curve(segmento,100);
        else
            inicio = index(1);
            fim = index(length(index))+1;
            ind = 0;
            A = length(segmento(1:inicio,:));
            B = length(segmento(fim:end,:));
            for k = 1:inicio
                prop = bm_length_of_curve(segmento,k)/...
                    (bm_length_of_curve(segmento,A)+ bm_length_of_curve(segmento(fim:end,:),B));
                if prop >= 0.5
                    ind = k;
                    break
                end
            end
            if ind == 0
                ind = floor(4*inicio/5);
            end
            midpoint = segmento(ind,:);
            equi_coord = [init_pos;midpoint;final_pos];
            handles.new_handle_order{end+1} = bm_handle(handles.num_of_handles+new_bones,equi_coord,mesh_coord,'Curved',0);
            handles.new_handle_order{end}.proportion = bm_length_of_curve(segmento,ind)/...
                (bm_length_of_curve(segmento,A)+ bm_length_of_curve(segmento(fim:end,:),B));
        end
        
    end
end

%% Here we take care of nearby endpoints
handles.num_of_handles = length(handles.new_handle_order);
handles.coupled_endpoints = cell(handles.num_of_handles,1);
handles.handles_visibility = ones(handles.num_of_handles,3);
handles.bone_intersections = zeros(handles.num_of_handles,handles.num_of_handles);
handles.intersection_points = [];
handles.point_coord_clusters = [];
handles.point_index_clusters = [];

for i  = 1:handles.num_of_handles
    if isequal(handles.new_handle_order{i}.type,'Curved')
        P1 = handles.new_handle_order{i}.position(1,:);
        P3 = handles.new_handle_order{i}.position(3,:);
        handles.point_coord_clusters = [handles.point_coord_clusters; P1 ; P3];
        handles.point_index_clusters = [handles.point_index_clusters; i 1; i 3];
    end
end
raio_approx = 2*pi/100;
[clustersCentroids,clustersGeoMedians,clustersXY, handles.clustersInd] = clusterXYpoints(handles.point_coord_clusters,raio_approx,1,'centroid','merge');

handles.num_of_clusters = numel(handles.clustersInd);
for i = 1:handles.num_of_clusters
    for j = handles.clustersInd{i}
        index = handles.point_index_clusters(j,:);
        handles.new_handle_order{index(1)}.position(index(2),:) = clustersCentroids(i,:);
        handles.new_handle_order{index(1)}.new_position = handles.new_handle_order{index(1)}.position;
        handles.new_handle_order{index(1)}.mesh_position = sphere2equi(handles.V(knnsearch(handles.V,equi2sphere(handles.new_handle_order{index(1)}.position)),:));
    end
end
for i = 1:handles.num_of_clusters
    if length(handles.clustersInd{i}) > 1
        for j = handles.clustersInd{i}(2:end)
            index = handles.point_index_clusters(j,:);
            handles.handles_visibility(index(1),index(2)) = 0;
        end
    end
end
for i = 1:length(handles.new_handle_order)
    if isequal(handles.new_handle_order{i}.type,'Curved')
        if isequal(handles.new_handle_order{i}.mesh_position(1,:),handles.new_handle_order{i}.mesh_position(3,:))
            handles.new_handle_order{i} = {};
        end
    end
end
handles.new_handle_order = handles.new_handle_order(~cellfun('isempty',handles.new_handle_order));
handles.num_of_handles = length(handles.new_handle_order);
handles.coupled_endpoints = cell(handles.num_of_handles,1);
handles.handles_visibility = ones(handles.num_of_handles,3);
handles.bone_intersections = zeros(handles.num_of_handles,handles.num_of_handles);
handles.intersection_points = [];
handles.point_coord_clusters = [];
handles.point_index_clusters = [];

for i  = 1:handles.num_of_handles
    if isequal(handles.new_handle_order{i}.type,'Curved')
        P1 = handles.new_handle_order{i}.position(1,:);
        P3 = handles.new_handle_order{i}.position(3,:);
        handles.point_coord_clusters = [handles.point_coord_clusters; P1 ; P3];
        handles.point_index_clusters = [handles.point_index_clusters; i 1; i 3];
    end
end
raio_approx = 2*pi/100;
[clustersCentroids,clustersGeoMedians,clustersXY, handles.clustersInd] = clusterXYpoints(handles.point_coord_clusters,raio_approx,1,'centroid','merge');

handles.num_of_clusters = numel(handles.clustersInd);
for i = 1:handles.num_of_clusters
    for j = handles.clustersInd{i}
        index = handles.point_index_clusters(j,:);
        handles.new_handle_order{index(1)}.position(index(2),:) = clustersCentroids(i,:);
        handles.new_handle_order{index(1)}.new_position = handles.new_handle_order{index(1)}.position;
        handles.new_handle_order{index(1)}.mesh_position = sphere2equi(handles.V(knnsearch(handles.V,equi2sphere(handles.new_handle_order{index(1)}.position)),:));
    end
end
for i = 1:handles.num_of_clusters
    if length(handles.clustersInd{i}) > 1
        for j = handles.clustersInd{i}(2:end)
            index = handles.point_index_clusters(j,:);
            handles.handles_visibility(index(1),index(2)) = 0;
        end
    end
end

point_indices = [];
curved_bone_indices = [];
positions = [];


for i = 1:handles.num_of_handles
    positions = cat(1,positions,handles.new_handle_order{i}.mesh_position);
end
for i = 1:handles.num_of_handles
    switch handles.new_handle_order{i}.type
        case 'Curved'
            indice1 = knnsearch(positions,handles.new_handle_order{i}.mesh_position(1,:));
            indice2 = knnsearch(positions,handles.new_handle_order{i}.mesh_position(3,:));
            curved_bone_indices = cat(1,curved_bone_indices,[indice1 indice2]);
        case 'Point'
            indice = knnsearch(positions,handles.new_handle_order{i}.mesh_position(1,:));
            point_indices = cat(1,point_indices,indice);
        case 'Closed_Cage'
            indice = knnsearch(positions,handles.new_handle_order{i}.mesh_position(1,:));
            point_indices = cat(1,point_indices,indice);
    end
end

%% Here we take care of conform regions
is_conform = zeros(size(handles.V,1),1);
if ~isempty(handles.region_positions)
    size_regions = size(handles.region_positions,1);
    equi_positions_regions = zeros(size_regions,8);
    
    equi_positions_regions(:,1) = handles.region_positions(:,1);
    equi_positions_regions(:,2) = handles.region_positions(:,2);
    
    equi_positions_regions(:,3) = handles.region_positions(:,1);
    equi_positions_regions(:,4) = handles.region_positions(:,2)+handles.region_positions(:,4);
    
    equi_positions_regions(:,5) = handles.region_positions(:,1)+handles.region_positions(:,3);
    equi_positions_regions(:,6) = handles.region_positions(:,2)+handles.region_positions(:,4);
    
    equi_positions_regions(:,7) = handles.region_positions(:,1)+handles.region_positions(:,3);
    equi_positions_regions(:,8) = handles.region_positions(:,2);
    
    for index_region = 1:2:7
        calc_pos = pixel2equi([equi_positions_regions(:,index_region),equi_positions_regions(:,index_region+1)],...
            handles.imagesize);
        equi_positions_regions(:,index_region) = calc_pos(:,1);
        equi_positions_regions(:,index_region+1) = calc_pos(:,2);
    end
    is_inside = zeros(size(handles.equirec_points,1),1);
    for other_index = 1:size_regions
        inside_polygon = inpolygon(handles.equirec_points(:,1),handles.equirec_points(:,2),...
            equi_positions_regions(other_index,1:2:7),equi_positions_regions(other_index,2:2:8));
        is_inside = max(is_inside,inside_polygon);
    end
    is_inside = [is_inside,is_inside];
    is_inside = is_inside.*handles.equirec_points;
    is_inside(is_inside(:,1)==0,:) = [];
    indices_equi = knnsearch(handles.V,equi2sphere(is_inside));
    is_conform(indices_equi) = 1;
end


% tic
% 
% [b,bc] = bm_boundary_conditions(handles.V,handles.F,equi2sphere(positions),point_indices,...
%     [],[],curved_bone_indices,finally_cage_indices);
% if ~isempty(is_conform)
%     W = biharmonic_bounded(handles.V,handles.F,b,bc,'POU',false,...
%         'ShapePreserving',is_conform);
% else
%     W = biharmonic_bounded(handles.V,handles.F,b,bc,'POU',false);
% end
% W = W(handles.indices,:);
% nx = handles.mesh_quality + 52;
% x = linspace(-pi+1e-6,pi-1e-6,nx);
% y = linspace(pi/2-1e-6,-pi/2+1e-6,nx/2);
% [X,Y]=meshgrid(x,y);
% 
% handles.non_expanded = []; %from irregular to regular mesh
% for i  = 1:handles.num_of_handles
%     handles.non_expanded(:,:,i) = griddata(handles.equirec_points(:,1),handles.equirec_points(:,2),W(:,i),X,Y);
% end
% nx = handles.imagesize(2);
% ny = handles.imagesize(1);
% x = linspace(-pi+1e-6,pi-1e-6,nx);
% y = linspace(pi/2-1e-6,-pi/2+1e-6,ny);
% [U,V]=meshgrid(x,y);
% handles.W_expanded = []; %expand to image size
% for i  = 1:handles.num_of_handles
%     handles.W_expanded(:,:,i) = qinterp2(X,Y,handles.non_expanded(:,:,i),U,V);
% end
% 
% time_string = toc;
% msg = strcat('Total computation time: ', num2str(time_string),' seconds');
% h = msgbox(msg,'Total time');
% uiwait(h);
transformation_interface(handles);

closereq;


% --- Executes on button press in load.
function load_Callback(hObject, eventdata, handles)
% hObject    handle to load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uigetfile({'*.jpg'},'File Selector');
inputfilename = strcat(pathname,filename);
handles.inputfilename = inputfilename;
handles.matlabImage = imread(handles.inputfilename);
handles.imagesize = size(handles.matlabImage);
axes(handles.image)
handles.img_handle = image(handles.matlabImage);
axis off

handles.num_of_handles = 0;
handles.all_handles = {};
handles.handles_plot = {};
handles.new_handle_order = {};

set(handles.img_handle,'ButtonDownFcn',{@image_ButtonDownFcn,handles});
set(handles.image,'ButtonDownFcn',{@image_ButtonDownFcn,handles});
guidata(hObject,handles)



% --- Executes when figure1 is resized.
function figure1_SizeChangedFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% LimitFigSize(handles.figure1, 'min', [840, 360])


window_values = get(handles.figure1,'Position');

%  image
Pos_image = get(handles.image,'Position');
Pos_image(3) = min([2*(window_values(4)-140) window_values(3)-40]);
Pos_image(4) = round(Pos_image(3)/2);
Pos_image(2) = round((window_values(4)-120)/2 - Pos_image(4)/2);
Pos_image(1) = (window_values(3)-Pos_image(3))/2;
set(handles.image,'Position',Pos_image);

% button_group
Pos_button = get(handles.uibuttongroup1,'Position');
Pos_button(1) = round(window_values(3)/2-320);
Pos_button(2) = round(window_values(4) - 120);
set(handles.uibuttongroup1,'Position',Pos_button);
% delete button
Pos_delete = get(handles.delete_handle,'Position');
Pos_delete(1) = round(window_values(3)/2 +20);
Pos_delete(2) = round(window_values(4) - 120);
set(handles.delete_handle,'Position',Pos_delete);

% next
Pos_next = get(handles.next,'Position');
Pos_next(1) = round(window_values(3)/2) + 240;
Pos_next(2) = round(window_values(4) - 60);
set(handles.next,'Position',Pos_next);

%delete region button
pos_region = get(handles.delete_region,'Position');
pos_region(1) = round(window_values(3)/2) + 200;
pos_region(2) = round(window_values(4) - 120);
set(handles.delete_region,'Position',pos_region);

%load
Pos_load = get(handles.load,'Position');
Pos_load(1) = round(window_values(3)/2 +20);
Pos_load(2) = round(window_values(4) - 60);
set(handles.load,'Position',Pos_load);



% --- Executes on button press in delete_region.
function delete_region_Callback(hObject, eventdata, handles)
% hObject    handle to delete_region (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isempty(handles.all_regions)
    delete(handles.all_regions{end})
    handles.all_regions(end) = [];
    handles.region_positions = handles.region_positions(1:end-1,:);
end
guidata(hObject,handles);


% --- Executes on button press in conformal_region.
function conformal_region_Callback(hObject, eventdata, handles)
% hObject    handle to conformal_region (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of conformal_region
