function varargout = new_interface(varargin)
% NEW_INTERFACE MATLAB code for new_interface.fig
%      NEW_INTERFACE, by itself, creates a new NEW_INTERFACE or raises the existing
%      singleton*.
%
%      H = NEW_INTERFACE returns the handle to a new NEW_INTERFACE or the handle to
%      the existing singleton*.
%
%      NEW_INTERFACE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in NEW_INTERFACE.M with the given input arguments.
%
%      NEW_INTERFACE('Property','Value',...) creates a new NEW_INTERFACE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before new_interface_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to new_interface_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help new_interface

% Last Modified by GUIDE v2.5 08-May-2018 20:51:27

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @new_interface_OpeningFcn, ...
                   'gui_OutputFcn',  @new_interface_OutputFcn, ...
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


% --- Executes just before new_interface is made visible.
function new_interface_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to new_interface (see VARARGIN)

% Choose default command line output for new_interface
handles.output = hObject;
handles.output = hObject;
handles.inputfilename = 'grabodrertorv.jpg';
handles.matlabImage = imread(handles.inputfilename);
handles.imagesize = size(handles.matlabImage);
morehandles = varargin{1};
handles.sphere_points = morehandles.sphere_points;
handles.equi_points = morehandles.equi_points;
handles.num_of_handles = morehandles.num_of_handles;
handles.W_expanded = morehandles.W_expanded;
handles.quality = morehandles.quality;
handles.F = morehandles.F;
handles.F2 = morehandles.F2;
handles.hd = {};
handles.all_handles = morehandles.all_handles;
handles.handles_plot = morehandles.handles_plot;
handles.handles_plot = {};
handles.current_handle = [];
% handles.mouse = vrspacemouse()
axes(handles.image)
handles.img_handle = image(handles.matlabImage);
hold on

for i= 1:handles.num_of_handles
    pix_pos = equi2pixels(handles.all_handles{i}.position,handles.imagesize);
    if isequal(handles.all_handles{i}.type,'Curved')
        pix_aux_pos = equi2pixels(handles.all_handles{i}.auxiliar_position,handles.imagesize);
        e = line('Parent', handles.image, 'XData', [pix_pos(1,1) pix_aux_pos(1) pix_pos(2,1)],...
             'YData', [pix_pos(1,2) pix_aux_pos(2) pix_pos(2,2)],'Color','b','LineStyle','--','LineWidth',2);
        f = plot(pix_pos(1,1),pix_pos(1,2),'ro','MarkerFaceColor','r','MarkerSize',4,'Parent',handles.image);
        g = plot(pix_pos(2,1),pix_pos(2,2),'ro','MarkerFaceColor','r','MarkerSize',4,'Parent',handles.image);
        h = plot(pix_aux_pos(1),pix_aux_pos(2),'ro','MarkerFaceColor','r','MarkerSize',4,'Parent',handles.image);
        handles.handles_plot{end+1} = {e,f,g,h};
    else
        pix_aux_pos = equi2pixels(handles.all_handles{i}.auxiliar_position,handles.imagesize);
        f = arrow([pix_pos(1) pix_pos(2)],[pix_aux_pos(1) pix_pos(2)],'Width',2,'Length',10,'FaceColor','b','EdgeColor','b');
        g = arrow([pix_pos(1) pix_pos(2)],[pix_pos(1)-(pix_aux_pos(1)-pix_pos(1)) pix_pos(2)],'Width',2,'Length',10,'FaceColor','g','EdgeColor','g');
        h = impoint(gca,pix_pos(1), pix_pos(2));
        handles.handles_plot{end+1} = {h,f,g};
        mycb1 = @(pos) Position_of_point(h,pos,handles,1);
        addNewPositionCallback(h,mycb1);
%         set(h,'ButtonDownFcn',{@release,handles});
        setColor(h,'red');
        
   end
end
axis off
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes new_interface wait for user response (see UIRESUME)
% uiwait(handles.figure1);
function Position_of_point(h,pos,handles,indice)

if numel(handles.handles_plot{1}) > 4
    delete(handles.handles_plot{1}{3})
    delete(handles.handles_plot{1}{2})
    handles.handles_plot{1} = handles.handles_plot{1}{1};
else
    handles.handles_plot{1}{2} = arrow([pos(1) pos(2)],[pos(1)+60 pos(2)],'Width',2,'Length',10,'FaceColor','b','EdgeColor','b');
	handles.handles_plot{1}{3} = arrow([pos(1) pos(2)],[pos(1)-(60) pos(2)],'Width',2,'Length',10,'FaceColor','g','EdgeColor','g');
end
%         mycb1 = @(pos) Position_of_point(h,pos,handles,1);
%         addNewPositionCallback(h,mycb1);

getPosition(handles.handles_plot{1}{1})

% function release(hObject,eventdata,handles)
% if numel(handles.handles_plot{1}) > 1
%     delete(handles.handles_plot{1}{3})
%     delete(handles.handles_plot{1}{2})
% end

% --- Outputs from this function are returned to the command line.
function varargout = new_interface_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
