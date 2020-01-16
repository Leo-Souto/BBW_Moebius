function varargout = plot_interface(varargin)
% PLOT_INTERFACE MATLAB code for plot_interface.fig
%      PLOT_INTERFACE, by itself, creates a new PLOT_INTERFACE or raises the existing
%      singleton*.
%
%      H = PLOT_INTERFACE returns the handle to a new PLOT_INTERFACE or the handle to
%      the existing singleton*.
%
%      PLOT_INTERFACE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PLOT_INTERFACE.M with the given input arguments.
%
%      PLOT_INTERFACE('Property','Value',...) creates a new PLOT_INTERFACE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before plot_interface_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to plot_interface_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help plot_interface

% Last Modified by GUIDE v2.5 03-Jun-2018 14:07:23

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @plot_interface_OpeningFcn, ...
                   'gui_OutputFcn',  @plot_interface_OutputFcn, ...
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


% --- Executes just before plot_interface is made visible.
function plot_interface_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to plot_interface (see VARARGIN)

% Choose default command line output for plot_interface
handles.output = hObject;
cor = load('custom_color.mat');
handles.color = cor.cor;
%Recuperando os handles da janela anterior
morehandles = varargin{1};
handles.sphere_points = morehandles.sphere_points;
handles.equi_points = morehandles.equirec_points;
handles.num_of_handles = morehandles.num_of_handles;
handles.W_expanded = morehandles.W_expanded;
handles.quality = morehandles.mesh_quality;
handles.inputfilename = morehandles.inputfilename;
handles.F = morehandles.F;
handles.F2 = morehandles.F2;
handles.all_cages = morehandles.new_cages;
handles.non_expanded = morehandles.non_expanded;
handles.bone_intersections = morehandles.bone_intersections;
handles.intersection_points = morehandles.intersection_points;
handles.coupled_endpoints = morehandles.coupled_endpoints;
handles.handles_visibility = morehandles.handles_visibility;
handles.point_coord_clusters = morehandles.point_coord_clusters;
handles.point_index_clusters = morehandles.point_index_clusters;
handles.clustersInd = morehandles.clustersInd;
handles.num_of_clusters = morehandles.num_of_clusters;
if handles.num_of_handles > 0
    handles.all_handles = morehandles.new_handle_order;
    handles.handles_plot = morehandles.handles_plot;
    c_string = {};
    for i = 1:handles.num_of_handles
    c_string{end+1} = num2str(i);
    end
    set(handles.handles_menu,'String',c_string)
    nx = handles.quality + 52;
    x = linspace(-pi,pi,nx);  
    y = linspace(pi/2,-pi/2,nx/2);  
    [handles.X, handles.Y]=meshgrid(x,y);
    A = 1;
axes(handles.plotter)

%     axis equal
%     hold on
%     colormap(handles.color)
    pcolor(handles.X,handles.Y,handles.non_expanded(:,:,A))
    colormap(handles.color)
    shading interp
    hold on
    contour(handles.X,handles.Y,handles.non_expanded(:,:,A),'k')
    view(0,90);
    axis off
    rotate3d off
    hold off
else
    warning('No handle available')
end
set(handles.v_number,'String',num2str(size(handles.equi_points,1)))
set(handles.f_number,'String',num2str(size(handles.F2,1)))

% Inicia o ponto
handles.point = [0 0];
handles.point2 = [0 0];
set(handles.figure1, 'Position', get(0, 'Screensize'));

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes plot_interface wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = plot_interface_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in handles_menu.
function handles_menu_Callback(hObject, eventdata, handles)
% hObject    handle to handles_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
A = get(hObject,'Value');
if ~isnan(A)

axes(handles.plotter)
    pcolor(handles.X,handles.Y,handles.non_expanded(:,:,A))
    colormap(handles.color)
    shading interp
    hold on
    contour(handles.X,handles.Y,handles.non_expanded(:,:,A),'k')
    view(0,90);
    axis off
    rotate3d off
    hold off
end
guidata(hObject,handles)
% Hints: contents = cellstr(get(hObject,'String')) returns handles_menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from handles_menu


% --- Executes during object creation, after setting all properties.
function handles_menu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to handles_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in next.
function next_Callback(hObject, eventdata, handles)
% hObject    handle to next (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
transformation_interface(handles);
% new_interface(handles);
closereq;


% --- Executes when figure1 is resized.
function figure1_SizeChangedFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
window_values = get(handles.figure1,'Position');

% plotter
Pos_plotter = get(handles.plotter,'Position');
Pos_plotter(2) = 40;
Pos_plotter(3) = min([2*(window_values(4)-120) window_values(3)-80]);
Pos_plotter(4) = round(Pos_plotter(3)/2);
Pos_plotter(1) = (window_values(3)-Pos_plotter(3))/2;
set(handles.plotter,'Position',Pos_plotter);


% vertices_text
Pos_vertices_text = get(handles.vertices_text,'Position');
Pos_vertices_text(1) = window_values(3)/2 - 290;
Pos_vertices_text(2) = window_values(4) - 40;
set(handles.vertices_text,'Position',Pos_vertices_text);
%v_number
Pos_v_number = get(handles.v_number,'Position');
Pos_v_number(1) = window_values(3)/2 - 290;
Pos_v_number(2) = window_values(4) - 60;
set(handles.v_number,'Position',Pos_v_number);

% faces_text
Pos_faces_text = get(handles.faces_text,'Position');
Pos_faces_text(1) = window_values(3)/2 - 90;
Pos_faces_text(2) = window_values(4) - 40;
set(handles.faces_text,'Position',Pos_faces_text);
%f_number
Pos_f_number = get(handles.f_number,'Position');
Pos_f_number(1) = window_values(3)/2 - 90;
Pos_f_number(2) = window_values(4) - 60;
set(handles.f_number,'Position',Pos_f_number);

%select
Pos_select = get(handles.select,'Position');
Pos_select(1) = window_values(3)/2 + 115;
Pos_select(2) = window_values(4) - 40;
set(handles.select,'Position',Pos_select);
%handles_menu
Pos_handles = get(handles.handles_menu,'Position');
Pos_handles(1) = window_values(3)/2 + 140;
Pos_handles(2) = window_values(4) - 60;
set(handles.handles_menu,'Position',Pos_handles);

% next button
pos_next = get(handles.next,'Position');
pos_next(1) = window_values(3)/2 + 220;
pos_next(2) = window_values(4)-60;
set(handles.next,'Position',pos_next);

