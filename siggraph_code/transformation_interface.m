function varargout = transformation_interface(varargin)
% TRANSFORMATION_INTERFACE MATLAB code for transformation_interface.fig
%      TRANSFORMATION_INTERFACE, by itself, creates a new TRANSFORMATION_INTERFACE or raises the existing
%      singleton*.
%
%      H = TRANSFORMATION_INTERFACE returns the handle to a new TRANSFORMATION_INTERFACE or the handle to
%      the existing singleton*.
%
%      TRANSFORMATION_INTERFACE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TRANSFORMATION_INTERFACE.M with the given input arguments.
%
%      TRANSFORMATION_INTERFACE('Property','Value',...) creates a new TRANSFORMATION_INTERFACE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before transformation_interface_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to transformation_interface_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help transformation_interface

% Last Modified by GUIDE v2.5 11-Jan-2019 13:39:10

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @transformation_interface_OpeningFcn, ...
    'gui_OutputFcn',  @transformation_interface_OutputFcn, ...
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


% --- Executes just before transformation_interface is made visible.
function transformation_interface_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to transformation_interface (see VARARGIN)

% Choose default command line output for transformation_interface
handles.output = hObject;
morehandles = varargin{1};
handles.inputfilename = morehandles.inputfilename;
handles.matlabImage = imread(handles.inputfilename);
handles.originalimage = imread(handles.inputfilename);
handles.imagesize = size(handles.matlabImage);
handles.V_equi = morehandles.V_equi;
handles.num_of_handles = morehandles.num_of_handles;
% handles.W = morehandles.W_expanded;
handles.all_handles = morehandles.all_handles;
handles.all_cages = morehandles.all_cages;
handles.bone_intersections = morehandles.bone_intersections;
handles.intersection_points = morehandles.intersection_points;
handles.coupled_endpoints = morehandles.coupled_endpoints;
handles.handles_visibility = morehandles.handles_visibility;
handles.point_coord_clusters = morehandles.point_coord_clusters;
handles.point_index_clusters = morehandles.point_index_clusters;
% handles.intersections_index = morehandles.intersections_index;
handles.clustersInd = morehandles.clustersInd;
handles.num_of_clusters = morehandles.num_of_clusters;
handles.V = morehandles.V;
handles.V_equi = morehandles.V_equi;
handles.F = morehandles.F;
handles.F_equi = morehandles.F_equi;
% handles.intersections_index = morehandles.intersections_index;
handles.handles_plot = {};
handles.current_handle = [];
handles.T_coefs = [];
handles.cross_end = [];
handles.plot_pos = {};
handles.cage_weights = {};
handles.cage_plot_pos = {};
handles.cage_cross_end = [];
handles.cages_plot = {};


handles.rectify_atual = zeros(handles.num_of_handles,1);
axes(handles.image)
handles.img_handle = image(handles.matlabImage);

hold on
x = linspace(-pi,pi,handles.imagesize(2));
y = linspace(-pi/2,pi/2,handles.imagesize(1));
[X,Y] = ndgrid(x,y);
for i = 1:length(handles.all_cages)
    for j = 1:length(handles.all_cages{i})-1
        index1 = handles.all_cages{i}(j);
        index2 = handles.all_cages{i}(j+1);
        interp_1 = griddedInterpolant(X,Y,(flipud(handles.W(:,:,index1)))');
        interp_2 = griddedInterpolant(X,Y,(flipud(handles.W(:,:,index2)))');
        sphere_coord = equi2sphere(handles.all_handles{index1}.position(1,:));
        sphere_coord2 = equi2sphere(handles.all_handles{index2}.position(1,:));
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
        weight1 = interp_1(equi_coord_par);
        weight2 = interp_2(equi_coord_par);
        handles.cage_weights{i}{j} = [weight1 weight2];
        deriv = diff(equi_coord_par(:,1));
        index = find(abs(deriv) > 0.1);
        index = sort(index);
        equi_coord_par_pix = equi2pixels(sphere2equi(parametric),handles.imagesize);
        handles.cage_plot_pos{i}{j} = sphere2equi(parametric);
        if ~isempty(index)
            inicio = index(1);
            fim = index(length(index))+1;
            X1 = equi_coord_par_pix(1:inicio,1);
            Y1 = equi_coord_par_pix(1:inicio,2);
            X2 = equi_coord_par_pix(fim:100,1);
            Y2 = equi_coord_par_pix(fim:100,2);
            f = plot(X1,Y1,'Color',[1 1 1],'Parent', handles.image);
            f2 = plot(X2,Y2,'Color',[1 1 1],'Parent', handles.image);
            handles.cages_plot{end+1} = f;
            handles.cages_plot{end+1} = f2;
            handles.cage_cross_end = [handles.cross_end i j];
        else
            f = plot(equi_coord_par_pix(:,1),equi_coord_par_pix(:,2),'Color',[1 1 1],'Parent', handles.image );
            handles.cages_plot{end+1} = f;
        end
    end
end


for i= 1:handles.num_of_handles
    pix_pos = equi2pixels(handles.all_handles{i}.position,handles.imagesize);
    equi_pos = handles.all_handles{i}.position;
    if isequal(handles.all_handles{i}.type,'Curved')
        sphere_coord = equi2sphere(equi_pos(1,:));
        sphere_coord2 = equi2sphere(equi_pos(3,:));
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
        handles.plot_pos{i} = equi_coord_par;
        deriv = diff(equi_coord_par(:,1));
        index = find(abs(deriv) > 0.1);
        index = sort(index);
        equi_coord_par_pix = equi2pixels(sphere2equi(parametric),handles.imagesize);
        if ~isempty(index)
            inicio = index(1);
            fim = index(length(index))+1;
            X1 = equi_coord_par_pix(1:inicio,1);
            Y1 = equi_coord_par_pix(1:inicio,2);
            X2 = equi_coord_par_pix(fim:100,1);
            Y2 = equi_coord_par_pix(fim:100,2);
            l1 = plot(X1,Y1,'Color',[0.79 0 0.125],'Parent', handles.image);
            l2 = plot(X2,Y2,'Color',[0.79 0 0.125],'Parent', handles.image);
            f =  plot(pix_pos(1,1),pix_pos(1,2),'o','MarkerFaceColor',[0.79 0 0.125],'MarkerEdgeColor',[0.25 0.25 0.25],'MarkerSize',5,'Parent',handles.image);
            g =  plot(pix_pos(2,1),pix_pos(2,2),'o','MarkerFaceColor',[0.79 0 0.125],'MarkerEdgeColor',[0.25 0.25 0.25],'MarkerSize',5,'Parent',handles.image);
            h =  plot(pix_pos(3,1),pix_pos(3,2),'o','MarkerFaceColor',[0.79 0 0.125],'MarkerEdgeColor',[0.25 0.25 0.25],'MarkerSize',5,'Parent',handles.image);
            handles.handles_plot{end+1} = {f,g,h,l1,l2};
        else
            l = plot(equi_coord_par_pix(:,1),equi_coord_par_pix(:,2),'Color',[0.79 0 0.125]);
            f =  plot(pix_pos(1,1),pix_pos(1,2),'o','MarkerFaceColor',[0.79 0 0.125],'MarkerEdgeColor',[0.25 0.25 0.25],'MarkerSize',5,'Parent',handles.image);
            g =  plot(pix_pos(2,1),pix_pos(2,2),'o','MarkerFaceColor',[0.79 0 0.125],'MarkerEdgeColor',[0.25 0.25 0.25],'MarkerSize',5,'Parent',handles.image);
            h =  plot(pix_pos(3,1),pix_pos(3,2),'o','MarkerFaceColor',[0.79 0 0.125],'MarkerEdgeColor',[0.25 0.25 0.25],'MarkerSize',5,'Parent',handles.image);
            handles.handles_plot{end+1} = {f,g,h,l};
        end
    elseif isequal(handles.all_handles{i}.type,'Point')
        f = arrow([pix_pos(1,1) pix_pos(1,2)],[pix_pos(2,1) pix_pos(2,2)],'Width',2,'Length',10,'FaceColor',[0.79 0 0.125],'EdgeColor',[0.25 0.25 0.25]);
        g = arrow([pix_pos(1,1) pix_pos(1,2)],[pix_pos(3,1) pix_pos(3,2)],'Width',2,'Length',10,'FaceColor',[0.79 0 0.125],'EdgeColor',[0.25 0.25 0.25]);
        h = plot(pix_pos(1,1),pix_pos(1,2),'o','MarkerFaceColor',[0.79 0 0.125],'MarkerEdgeColor',[0.25 0.25 0.25],'MarkerSize',5,'Parent',handles.image);
        handles.handles_plot{end+1} = {h,f,g};
    elseif isequal(handles.all_handles{i}.type,'Closed_Cage')
        f = plot(pix_pos(1,1),pix_pos(1,2),'o','MarkerFaceColor',[0.79 0 0.125],'MarkerEdgeColor',[0.25 0.25 0.25],'MarkerSize',5,'Parent',handles.image);
        handles.handles_plot{end+1} = f;
    end
    handles.T_coefs(i,:) = [1 0 0 1];
end
for i = 1:handles.num_of_handles
    if isequal(handles.all_handles{i}.type,'Point')
        set(handles.handles_plot{i}{1},'UserData',[i 1]);
        set(handles.handles_plot{i}{2},'UserData',[i 2]);
        set(handles.handles_plot{i}{3},'UserData',[i 3]);
        set(handles.handles_plot{i}{1},'ButtonDownFcn',{@selectpoint,handles})
        set(handles.handles_plot{i}{2},'ButtonDownFcn',{@selectpoint,handles})
        set(handles.handles_plot{i}{3},'ButtonDownFcn',{@selectpoint,handles})
    elseif isequal(handles.all_handles{i}.type,'Curved')
        set(handles.handles_plot{i}{1},'UserData',[i 1]);
        set(handles.handles_plot{i}{2},'UserData',[i 2]);
        set(handles.handles_plot{i}{3},'UserData',[i 3]);
        set(handles.handles_plot{i}{1},'ButtonDownFcn',{@selectpoint,handles})
        set(handles.handles_plot{i}{2},'ButtonDownFcn',{@selectpoint,handles})
        set(handles.handles_plot{i}{3},'ButtonDownFcn',{@selectpoint,handles})
    elseif isequal(handles.all_handles{i}.type,'Closed_Cage')
        set(handles.handles_plot{i},'UserData',[i 1]);
        set(handles.handles_plot{i},'ButtonDownFcn',{@selectpoint,handles})
    end
end
for i = 1:handles.num_of_handles
    if isequal(handles.all_handles{i}.type,'Curved')
        
        uistack(handles.handles_plot{i}{1},'top')
        uistack(handles.handles_plot{i}{2},'top')
        uistack(handles.handles_plot{i}{3},'top')
        if ~handles.handles_visibility(i,1)
            set(handles.handles_plot{i}{1},'Visible','off');
        end
        if ~handles.handles_visibility(i,3)
            set(handles.handles_plot{i}{3},'Visible','off');
        end
    end
end
axis off
hold off

set(handles.figure1, 'Position', get(0, 'Screensize'));

% Update handles structure
setappdata(0,'Handles',handles);
set(handles.img_handle, 'ButtonDownFcn', {@image_ButtonDownFcn,handles});
guidata(hObject, handles);

% UIWAIT makes transformation_interface wait for user response (see UIRESUME)
% uiwait(handles.figure1);
function image_ButtonDownFcn(hObject,eventdata,handles)
handles = getappdata(0,'Handles');
pt = get(get(hObject, 'Parent'), 'CurrentPoint');
x = pt(1,1);
y = pt(1,2);
i = handles.current_handle(1);
iii = handles.current_handle(2);

P_new = pixel2equi([x y],handles.imagesize);

if isequal(handles.all_handles{i}.type,'Point')
    P1 = handles.all_handles{i}.position(1,:);
    P2 = handles.all_handles{i}.position(2,:);
    P3 = handles.all_handles{i}.position(3,:);
    switch iii
        case 1
            xlength = handles.all_handles{i}.new_position(2,1) - handles.all_handles{i}.new_position(1,1);
            ylength = handles.all_handles{i}.new_position(2,2) - handles.all_handles{i}.new_position(1,2);
            handles.all_handles{i}.new_position(1,:) = P_new;
            P1_n = P_new;
            P2_n = [P1_n(1)+xlength P1_n(2)+ylength];
            P3_n = [P1_n(1)-xlength P1_n(2)-ylength];
        case 2
            xlength = P_new(1) - handles.all_handles{i}.new_position(1,1);
            ylength = P_new(2) - handles.all_handles{i}.new_position(1,2);
            P1_n = handles.all_handles{i}.new_position(1,:);
            P2_n = [P1_n(1)+xlength P1_n(2)+ylength];
            P3_n = [P1_n(1)-xlength P1_n(2)-ylength];
        case 3
            xlength = P_new(1) - handles.all_handles{i}.new_position(1,1);
            ylength = P_new(2) - handles.all_handles{i}.new_position(1,2);
            P1_n = handles.all_handles{i}.new_position(1,:);
            P2_n = [P1_n(1)-xlength P1_n(2)-ylength];
            P3_n = [P1_n(1)+xlength P1_n(2)+ylength];
    end
    handles.all_handles{i}.new_position = [P1_n;P2_n;P3_n];
    handles.T_coefs(i,:) = Calculate_Transformation([P1;P2;P3],[P1_n;P2_n;P3_n]);
elseif isequal(handles.all_handles{i}.type,'Curved')
    P1 = handles.all_handles{i}.position(1,:);
    P2 = handles.all_handles{i}.position(2,:);
    P3 = handles.all_handles{i}.position(3,:);
    handles.all_handles{i}.new_position(iii,:) = P_new;
    P1_n = handles.all_handles{i}.new_position(1,:);
    P3_n = handles.all_handles{i}.new_position(3,:);
    switch get(handles.rectify,'Value')
        case 1
            prop = handles.all_handles{i}.proportion;
            equi_coord_par = apply_transf_plot(handles.plot_pos{i},'Curved',handles.T_coefs(i,:),1);
            deriv = diff(equi_coord_par(:,1));
            index = find(abs(deriv) > 0.1);
            index = sort(index);
            if ~isempty(index)
                if (P1_n(1) - P3_n(1)) > P1_n(1)
                    P2_n = [P1_n(1)*(1-prop) + (prop)*(P3_n(1)+2*pi) P1_n(2)*(1-prop)+(prop)*P3_n(2)] ;
                elseif (P1_n(1) - P3_n(1))< P1_n(1)
                    P2_n = [P1_n(1)*(1-prop) + prop*(P3_n(1)+2*pi) P1_n(2)*(1-prop)+prop*P3_n(2)] ;
                end
            else
                P2_n = P1_n*(1-prop) + prop*P3_n;
            end
            
        case 0
            if handles.rectify_atual(i) == 1
                prop = handles.all_handles{i}.proportion;
                equi_coord_par = apply_transf_plot(handles.plot_pos{i},'Curved',handles.T_coefs(i,:),1);
                deriv = diff(equi_coord_par(:,1));
                index = find(abs(deriv) > 0.1);
                index = sort(index);
                if ~isempty(index)
                    if (P1_n(1) - P3_n(1)) > P1_n(1)
                        P2_n = [P1_n(1)*(1-prop) + prop*(P3_n(1)+2*pi) P1_n(2)*(1-prop)+prop*P3_n(2)] ;
                    elseif (P1_n(1) - P3_n(1))< P1_n(1)
                        P2_n = [P1_n(1)*(1-prop) + prop*(P3_n(1)-2*pi) P1_n(2)*(1-prop)+prop*P3_n(2)] ;
                    end
                else
                    P2_n = P1_n*(1-prop) + prop*P3_n;
                end
            else
                P2_n = handles.all_handles{i}.new_position(2,:);
            end
    end
    handles.all_handles{i}.new_position = [P1_n;P2_n;P3_n];
    handles.T_coefs(i,:) = Calculate_Transformation([P1;P2;P3],[P1_n;P2_n;P3_n]);
    handles.all_handles{i}.new_mesh_position = apply_transf_plot(handles.all_handles{i}.mesh_position,'Curved',handles.T_coefs(i,:),1);
    if iii ~= 3
        index_row = find(ismember(handles.point_index_clusters,[i iii],'rows'));
        index_cluster = [];
        for j = 1:handles.num_of_clusters
            if find(handles.clustersInd{j} == index_row)
                index_cluster = j;
            end
        end
        
        if ~isempty(index_cluster)
            for j = handles.clustersInd{index_cluster}
                i_c_p = handles.point_index_clusters(j,:);
                handles.all_handles{i_c_p(1)}.new_position(i_c_p(2),:) = P_new;
                P1_n = handles.all_handles{i_c_p(1)}.new_position(1,:);
                P3_n = handles.all_handles{i_c_p(1)}.new_position(3,:);
                switch get(handles.rectify,'Value')
                    case 1
                        prop = handles.all_handles{i_c_p(1)}.proportion;
                        equi_coord_par = apply_transf_plot(handles.plot_pos{i_c_p(1)},'Curved',handles.T_coefs(i_c_p(1),:),1);
                        deriv = diff(equi_coord_par(:,1));
                        index = find(abs(deriv) > 0.1);
                        index = sort(index);
                        if ~isempty(index)
                            if (P1_n(1) - P3_n(1))> P1_n(1)
                                P2_n = [P1_n(1)*(1-prop) + prop*(P3_n(1)+2*pi) P1_n(2)*(1-prop)+prop*P3_n(2)] ;
                            elseif (P1_n(1) - P3_n(1))< P1_n(1)
                                P2_n = [P1_n(1)*(1-prop) + prop*(P3_n(1)+2*pi) P1_n(2)*(1-prop)+prop*P3_n(2)] ;
                            end
                        else
                            P2_n = P1_n*(1-prop) + prop*P3_n;
                        end
                    case 0
                        if handles.rectify_atual(i_c_p(1)) == 1
                            prop = handles.all_handles{i_c_p(1)}.proportion;
                            equi_coord_par = apply_transf_plot(handles.plot_pos{i_c_p(1)},'Curved',handles.T_coefs(i_c_p(1),:),1);
                            deriv = diff(equi_coord_par(:,1));
                            index = find(abs(deriv) > 0.1);
                            index = sort(index);
                            if ~isempty(index)
                                if (P1_n(1) - P3_n(1))> P1_n(1)
                                    P2_n = [P1_n(1)*(1-prop) + prop*(P3_n(1)+2*pi) P1_n(2)*(1-prop)+prop*P3_n(2)] ;
                                elseif (P1_n(1) - P3_n(1))< P1_n(1)
                                    P2_n = [P1_n(1)*(1-prop) + prop*(P3_n(1)-2*pi) P1_n(2)*(1-prop)+prop*P3_n(2)] ;
                                end
                            else
                                P2_n = P1_n*(1-prop) + prop*P3_n;
                            end
                        else
                            P2_n = handles.all_handles{i_c_p(1)}.new_position(2,:);
                        end
                end
                handles.all_handles{i_c_p(1)}.new_position = [P1_n;P2_n;P3_n];
                handles.T_coefs(i_c_p(1),:) = Calculate_Transformation(handles.all_handles{i_c_p(1)}.position,handles.all_handles{i_c_p(1)}.new_position);
            end
        end
    end
elseif isequal(handles.all_handles{i}.type,'Closed_Cage')
    P1 = handles.all_handles{i}.position(1,:);
    P2 = handles.all_handles{i}.position(2,:);
    P3 = handles.all_handles{i}.position(3,:);
    P1_n = P_new;
    P2_n = [P_new(1)+0.3 P_new(2)];
    P3_n = [P_new(1)-0.3 P_new(2)];
    handles.all_handles{i}.new_position = [P1_n;P2_n;P3_n];
    handles.T_coefs(i,:) = Calculate_Transformation([P1;P2;P3],[P1_n;P2_n;P3_n]);
end

handles.matlabImage = biharmonic_moebius_sphere(handles.originalimage,handles.V,handles.F,handles.all_handles,handles.T_coefs,...
    [handles.imagesize(1) handles.imagesize(2)]);
axes(handles.image)
cla
handles.img_handle = image(handles.matlabImage);
hold on
handles.handles_plot = {};
handles.cages_plot = {};
for k = 1:length(handles.all_cages)
    for j = 1:length(handles.all_cages{k})-1
        index1 = handles.all_cages{k}(j);
        index2 = handles.all_cages{k}(j+1);
        Weight = handles.cage_weights{k}{j};
        transf = [handles.T_coefs(index1,:);handles.T_coefs(index2,:)];
        position = handles.cage_plot_pos{k}{j};
        equi_coord_par = apply_transf_plot(position,'Closed_Cage',transf,Weight);
        deriv = diff(equi_coord_par(:,1));
        index = find(abs(deriv) > 0.1);
        index = sort(index);
        equi_coord_par_pix = equi2pixels(equi_coord_par,handles.imagesize);
        if ~isempty(index)
            inicio = index(1);
            fim = index(length(index))+1;
            X1 = equi_coord_par_pix(1:inicio,1);
            Y1 = equi_coord_par_pix(1:inicio,2);
            X2 = equi_coord_par_pix(fim:100,1);
            Y2 = equi_coord_par_pix(fim:100,2);
            f = plot(X1,Y1,'Color',[1 1 1],'Parent', handles.image);
            f2 = plot(X2,Y2,'Color',[1 1 1],'Parent', handles.image);
            handles.cages_plot{end+1} = f;
            handles.cages_plot{end+1} = f2;
        else
            f = plot(equi_coord_par_pix(:,1), equi_coord_par_pix(:,2),'Color',[1 1 1],'Parent', handles.image);
            handles.cages_plot{end+1} = f;
        end
    end
end
for j= 1:handles.num_of_handles
    pix_pos = equi2pixels(handles.all_handles{j}.new_position,handles.imagesize);
    if isequal(handles.all_handles{j}.type,'Curved')
        equi_coord_par = apply_transf_plot(handles.plot_pos{j},'Curved',handles.T_coefs(j,:),1);
        equi_coord_par_pix = equi2pixels(equi_coord_par,handles.imagesize);
        deriv = diff(equi_coord_par(:,1));
        index = find(abs(deriv) > 0.1);
        index = sort(index);
        if ~isempty(index)
            inicio = index(1);
            fim = index(length(index))+1;
            X1 = equi_coord_par_pix(1:inicio,1);
            Y1 = equi_coord_par_pix(1:inicio,2);
            X2 = equi_coord_par_pix(fim:100,1);
            Y2 = equi_coord_par_pix(fim:100,2);
            l1 = plot(X1,Y1,'Color',[0.79 0 0.125],'Parent', handles.image);
            l2 = plot(X2,Y2,'Color',[0.79 0 0.125],'Parent', handles.image);
            f =  plot(pix_pos(1,1),pix_pos(1,2),'o','MarkerFaceColor',[0.79 0 0.125],'MarkerEdgeColor',[0.25 0.25 0.25],'MarkerSize',5,'Parent',handles.image);
            g =  plot(pix_pos(2,1),pix_pos(2,2),'o','MarkerFaceColor',[0.79 0 0.125],'MarkerEdgeColor',[0.25 0.25 0.25],'MarkerSize',5,'Parent',handles.image);
            h =  plot(pix_pos(3,1),pix_pos(3,2),'o','MarkerFaceColor',[0.79 0 0.125],'MarkerEdgeColor',[0.25 0.25 0.25],'MarkerSize',5,'Parent',handles.image);
            handles.handles_plot{end+1} = {f,g,h,l1,l2};
        else
            l = plot(equi_coord_par_pix(:,1),equi_coord_par_pix(:,2),'Color',[0.79 0 0.125]);
            f =  plot(pix_pos(1,1),pix_pos(1,2),'o','MarkerFaceColor',[0.79 0 0.125],'MarkerEdgeColor',[0.25 0.25 0.25],'MarkerSize',5,'Parent',handles.image);
            g =  plot(pix_pos(2,1),pix_pos(2,2),'o','MarkerFaceColor',[0.79 0 0.125],'MarkerEdgeColor',[0.25 0.25 0.25],'MarkerSize',5,'Parent',handles.image);
            h =  plot(pix_pos(3,1),pix_pos(3,2),'o','MarkerFaceColor',[0.79 0 0.125],'MarkerEdgeColor',[0.25 0.25 0.25],'MarkerSize',5,'Parent',handles.image);
            handles.handles_plot{end+1} = {f,g,h,l};
        end
        if handles.rectify_atual(j) == 1
            set(handles.handles_plot{j}{2},'Visible','off');
        end
        if ~handles.handles_visibility(j,1)
            set(handles.handles_plot{j}{1},'Visible','off');
        end
        if ~handles.handles_visibility(j,2)
            set(handles.handles_plot{j}{2},'Visible','off');
        end
        if ~handles.handles_visibility(j,3)
            set(handles.handles_plot{j}{3},'Visible','off');
        end
        set(handles.handles_plot{j}{1},'ButtonDownFcn',{@selectpoint,handles})
        set(handles.handles_plot{j}{2},'ButtonDownFcn',{@selectpoint,handles})
        set(handles.handles_plot{j}{3},'ButtonDownFcn',{@selectpoint,handles})
    elseif isequal(handles.all_handles{j}.type,'Point')
        f = arrow([pix_pos(1,1) pix_pos(1,2)],[pix_pos(2,1) pix_pos(2,2)],'Width',2,'Length',10,'FaceColor',[0.79 0 0.125],'EdgeColor',[0.25 0.25 0.25]);
        g = arrow([pix_pos(1,1) pix_pos(1,2)],[pix_pos(3,1) pix_pos(3,2)],'Width',2,'Length',10,'FaceColor',[0.79 0 0.125],'EdgeColor',[0.25 0.25 0.25]);
        h = plot(pix_pos(1,1),pix_pos(1,2),'o','MarkerFaceColor',[0.79 0 0.125],'MarkerEdgeColor',[0.25 0.25 0.25],'MarkerSize',5,'Parent',handles.image);
        handles.handles_plot{end+1} = {h,f,g};
        set(handles.handles_plot{j}{1},'UserData',[j 1]);
        set(handles.handles_plot{j}{2},'UserData',[j 2]);
        set(handles.handles_plot{j}{3},'UserData',[j 3]);
        set(handles.handles_plot{j}{1},'ButtonDownFcn',{@selectpoint,handles})
        set(handles.handles_plot{j}{2},'ButtonDownFcn',{@selectpoint,handles})
        set(handles.handles_plot{j}{3},'ButtonDownFcn',{@selectpoint,handles})
    elseif isequal(handles.all_handles{j}.type,'Closed_Cage')
        f = plot(pix_pos(1,1),pix_pos(1,2),'o','MarkerFaceColor',[0.79 0 0.125],'MarkerEdgeColor',[0.25 0.25 0.25],'MarkerSize',5,'Parent',handles.image);
        handles.handles_plot{j} = f;
        set(handles.handles_plot{j},'UserData',[j 1]);
        set(handles.handles_plot{j},'ButtonDownFcn',{@selectpoint,handles})
    end
end
axis off

if isequal(handles.all_handles{i}.type,'Point') && handles.current_handle(1,2) ~= 1
    set(handles.handles_plot{i}{iii},'FaceColor','y','EdgeColor',[0.25 0.25 0.25])
    set(handles.handles_plot{i}{iii}, 'ButtonDownFcn', {@image_ButtonDownFcn,handles});
elseif isequal(handles.all_handles{i}.type,'Closed_Cage')
    set(handles.handles_plot{i},'MarkerFaceColor','y','MarkerEdgeColor',[0.25 0.25 0.25])
else
    set(handles.handles_plot{i}{iii},'MarkerFaceColor','y','MarkerEdgeColor',[0.25 0.25 0.25])
end
hold off
setappdata(0,'Handles',handles);
hide_Callback(handles.hide,eventdata,handles);
Atualiza_handles();
setappdata(0,'Handles',handles);
% guidata(hObject,handles);

function new_pos = apply_transf_plot(position,type_flag,transf_coef,weights)
%% transforma handle
sphere_coord = equi2sphere(position);
%From sphere to stereographic plane
X_C =(2*sphere_coord(:,2))./(sphere_coord(:,1) + 1);
Y_C =(2*sphere_coord(:,3))./(sphere_coord(:,1) + 1);
%From stereographic plane to complex numbers
Complex = X_C + Y_C*1i;
Delta = zeros(size(sphere_coord));
if isequal(type_flag,'Curved')
    Complex2 = (transf_coef(1)*Complex + transf_coef(2))./(transf_coef(3)*Complex + transf_coef(4));
    X_C_n = real(Complex2);
    Y_C_n = imag(Complex2);
    
    r_polar = sqrt(X_C_n.^2 + Y_C_n.^2);
    theta_polar = atan2(X_C_n,Y_C_n);
    theta_polar = theta_polar + pi/2;
    X_C_n = r_polar.*cos(theta_polar);
    Y_C_n = r_polar.*sin(theta_polar);
    
    x_sphere2 = -(X_C_n.^2 + Y_C_n.^2 - 4)./(X_C_n.^2 + Y_C_n.^2 + 4);
    y_sphere2 = -(4*X_C_n)./(X_C_n.^2 + Y_C_n.^2 + 4);
    z_sphere2 = (4*Y_C_n)./(X_C_n.^2 + Y_C_n.^2 + 4);
    sphere_final = [x_sphere2  y_sphere2 z_sphere2];
    
elseif isequal(type_flag,'Closed_Cage')
    for i = [1 2]
        Complex2 = (transf_coef(i,1)*Complex + transf_coef(i,2))./(transf_coef(i,3)*Complex + transf_coef(i,4));
        %         Complex2 = exp(-pi/2*1i)*Complex2;
        X_C_n = real(Complex2);
        Y_C_n = imag(Complex2);
        
        r_polar = sqrt(X_C_n.^2 + Y_C_n.^2);
        theta_polar = atan2(X_C_n,Y_C_n);
        theta_polar = theta_polar + pi/2;
        X_C_n = r_polar.*cos(theta_polar);
        Y_C_n = r_polar.*sin(theta_polar);
        x_sphere2 = -(X_C_n.^2 + Y_C_n.^2 - 4)./(X_C_n.^2 + Y_C_n.^2 + 4);
        y_sphere2 = -(4*X_C_n)./(X_C_n.^2 + Y_C_n.^2 + 4);
        z_sphere2 = (4*Y_C_n)./(X_C_n.^2 + Y_C_n.^2 + 4);
        sphere2 = [x_sphere2  y_sphere2 z_sphere2];
        
        Delta(:,1) = Delta(:,1) - (weights(:,i)).*(sphere_coord(:,1) - sphere2(:,1));
        Delta(:,2) = Delta(:,2) - (weights(:,i)).*(sphere_coord(:,2) - sphere2(:,2));
        Delta(:,3) = Delta(:,3) - (weights(:,i)).*(sphere_coord(:,3) - sphere2(:,3));
    end
    sphere_final = sphere_coord + Delta;
end
new_pos = sphere2equi(sphere_final);

function selectpoint(hObject,eventdata,handles)
handles = getappdata(0,'Handles');
handles.current_handle = get(hObject,'UserData');
for i = 1:handles.num_of_handles
    if isequal(handles.all_handles{i}.type,'Point')
        set(handles.handles_plot{i}{1},'MarkerFaceColor',[0.79 0 0.125],'MarkerEdgeColor',[0.25 0.25 0.25])
        set(handles.handles_plot{i}{2},'FaceColor',[0.79 0 0.125],'EdgeColor',[0.25 0.25 0.25])
        set(handles.handles_plot{i}{3},'FaceColor',[0.79 0 0.125],'EdgeColor',[0.25 0.25 0.25])
    elseif isequal(handles.all_handles{i}.type,'Curved')
        set(handles.handles_plot{i}{1},'MarkerFaceColor',[0.79 0 0.125],'MarkerEdgeColor',[0.25 0.25 0.25])
        set(handles.handles_plot{i}{2},'MarkerFaceColor',[0.79 0 0.125],'MarkerEdgeColor',[0.25 0.25 0.25])
        set(handles.handles_plot{i}{3},'MarkerFaceColor',[0.79 0 0.125],'MarkerEdgeColor',[0.25 0.25 0.25])
    elseif isequal(handles.all_handles{i}.type,'Closed_Cage')
        set(handles.handles_plot{i},'MarkerFaceColor',[0.79 0 0.125],'MarkerEdgeColor',[0.25 0.25 0.25])
    end
end

setappdata(0,'Handles',handles);
Atualiza_handles()
handles = getappdata(0,'Handles');
i = handles.current_handle(1);
iii = handles.current_handle(2);
if handles.rectify_atual(i) == 1
    set(handles.rectify_current_line,'Value',1)
else
    set(handles.rectify_current_line,'Value',0)
end
if isequal(handles.all_handles{i}.type,'Point') && iii ~= 1
    set(handles.handles_plot{i}{iii},'FaceColor','y','EdgeColor',[0.25 0.25 0.25])
    set(hObject, 'ButtonDownFcn', {@image_ButtonDownFcn,handles});
else
    set(hObject,'MarkerFaceColor','y','MarkerEdgeColor',[0.25 0.25 0.25])
end
setappdata(0,'Handles',handles);


function Atualiza_handles()
handles = getappdata(0,'Handles');
for i = 1:handles.num_of_handles
    if isequal(handles.all_handles{i}.type,'Point')
        A = get(handles.handles_plot{i}{1},'UserData');
        B = get(handles.handles_plot{i}{2},'UserData');
        C = get(handles.handles_plot{i}{3},'UserData');
        if isempty(A)
            set(handles.handles_plot{i}{1},'UserData',[i 1]);
        end
        if isempty(B)
            set(handles.handles_plot{i}{2},'UserData',[i 2]);
        end
        if isempty(C)
            set(handles.handles_plot{i}{3},'UserData',[i 3]);
        end
        set(handles.handles_plot{i}{1},'ButtonDownFcn',{@selectpoint,handles})
        set(handles.handles_plot{i}{2},'ButtonDownFcn',{@selectpoint,handles})
        set(handles.handles_plot{i}{3},'ButtonDownFcn',{@selectpoint,handles})
    elseif isequal(handles.all_handles{i}.type,'Curved')
        A = get(handles.handles_plot{i}{1},'UserData');
        B = get(handles.handles_plot{i}{2},'UserData');
        C = get(handles.handles_plot{i}{3},'UserData');
        if isempty(A)
            set(handles.handles_plot{i}{1},'UserData',[i 1]);
        end
        if isempty(B)
            set(handles.handles_plot{i}{2},'UserData',[i 2]);
        end
        if isempty(C)
            set(handles.handles_plot{i}{3},'UserData',[i 3]);
        end
        set(handles.handles_plot{i}{1},'ButtonDownFcn',{@selectpoint,handles})
        set(handles.handles_plot{i}{2},'ButtonDownFcn',{@selectpoint,handles})
        set(handles.handles_plot{i}{3},'ButtonDownFcn',{@selectpoint,handles})
    elseif isequal(handles.all_handles{i}.type,'Closed_Cage')
        set(handles.handles_plot{i},'UserData',[i 1]);
        set(handles.handles_plot{i},'ButtonDownFcn',{@selectpoint,handles})
    end
end
set(handles.img_handle, 'ButtonDownFcn', {@image_ButtonDownFcn,handles});

setappdata(0,'Handles',handles);

function transf = Calculate_Transformation(inputs,outputs)
%Transform the input points in sphere points
input_coord = cat(2, inputs, outputs);
sphere_coord = zeros(3,6);
n = 1;
N = 3;
while n <= N
    sphere_coord(n,1) = cos(input_coord(n,2))*cos(input_coord(n,1));
    sphere_coord(n,2) = cos(input_coord(n,2))*sin(input_coord(n,1));
    sphere_coord(n,3) = sin(input_coord(n,2));
    sphere_coord(n,4) = cos(input_coord(n,4))*cos(input_coord(n,3));
    sphere_coord(n,5) = cos(input_coord(n,4))*sin(input_coord(n,3));
    sphere_coord(n,6) = sin(input_coord(n,4));
    n = n + 1;
end

%Projects the spheric input points onto stereografic plane
stereo_coord = zeros(3,4);
n = 1;
while n <= N
    stereo_coord(n,1) = 2*sphere_coord(n,2)/(sphere_coord(n,1) +1);
    stereo_coord(n,2) = 2*sphere_coord(n,3)/(sphere_coord(n,1) +1);
    stereo_coord(n,3) = 2*sphere_coord(n,5)/(sphere_coord(n,4) +1);
    stereo_coord(n,4) = 2*sphere_coord(n,6)/(sphere_coord(n,4) +1);
    n = n + 1;
end
%transform the stereographic input points onto complex numbers
complex_coord = zeros(3, 2);
n = 1;
while n <= N
    complex_coord(n,1) = stereo_coord(n,1) + stereo_coord(n,2)*1i;
    complex_coord(n,2) = stereo_coord(n,3) + stereo_coord(n,4)*1i;
    n = n + 1;
end
%Calculates the parameters of moebius transformation using a external
%function called crossmoebius
[a ,b, c, d] = crossmoebius(complex_coord(1,1), complex_coord(2,1), complex_coord(3,1),...
    complex_coord(1,2), complex_coord(2,2), complex_coord(3,2));
if abs(c) < 1e-10
    c = 1e-10 + 1e-10*1i;
end
transf = [a b c d];

% --- Executes on button press in undo.
function undo_Callback(hObject, eventdata, handles)
% hObject    handle to undo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = getappdata(0,'Handles');
if ~isempty(handles.current_handle)
    i = handles.current_handle(1);
    j = handles.current_handle(2);
    if isequal(handles.all_handles{i}.type,'Point')
        switch j
            case 1
                xlength = handles.all_handles{i}.new_position(2,1) - handles.all_handles{i}.new_position(1,1);
                ylength = handles.all_handles{i}.new_position(2,2) - handles.all_handles{i}.new_position(1,2);
                P1_n = handles.all_handles{i}.position(1,:);
                P2_n = [P1_n(1)+xlength P1_n(2)+ylength];
                P3_n = [P1_n(1)-xlength P1_n(2)-ylength];
            case 2
                xlength = 0.3;
                ylength = 0.0;
                P1_n = handles.all_handles{i}.new_position(1,:);
                P2_n = [P1_n(1)+xlength P1_n(2)+ylength];
                P3_n = [P1_n(1)-xlength P1_n(2)-ylength];
            case 3
                xlength = 0.3;
                ylength = 0.0;
                P1_n = handles.all_handles{i}.new_position(1,:);
                P2_n = [P1_n(1)+xlength P1_n(2)+ylength];
                P3_n = [P1_n(1)-xlength P1_n(2)-ylength];
        end
        P1 = handles.all_handles{i}.position(1,:);
        P2 = handles.all_handles{i}.position(2,:);
        P3 = handles.all_handles{i}.position(3,:);
        handles.all_handles{i}.new_position = [P1_n;P2_n;P3_n];
        handles.T_coefs(i,:) = Calculate_Transformation([P1;P2;P3],[P1_n;P2_n;P3_n]);
    elseif isequal(handles.all_handles{i}.type,'Curved')
        handles.all_handles{i}.new_position(j,:) = handles.all_handles{i}.position(j,:);
        P1_n = handles.all_handles{i}.new_position(1,:);
        P3_n = handles.all_handles{i}.new_position(3,:);
        switch get(handles.rectify,'Value')
            case 1
                prop = handles.all_handles{i}.proportion;
                equi_coord_par = apply_transf_plot(handles.plot_pos{i},'Curved',handles.T_coefs(i,:),1);
                deriv = diff(equi_coord_par(:,1));
                index = find(abs(deriv) > 0.1);
                index = sort(index);
                if ~isempty(index)
                    if (P1_n(1) - P3_n(1)) > P1_n(1)
                        P2_n = [P1_n(1)*(1-prop) + (prop)*(P3_n(1)+2*pi) P1_n(2)*(1-prop)+(prop)*P3_n(2)] ;
                    elseif (P1_n(1) - P3_n(1))< P1_n(1)
                        P2_n = [P1_n(1)*(1-prop) + prop*(P3_n(1)+2*pi) P1_n(2)*(1-prop)+prop*P3_n(2)] ;
                    end
                else
                    P2_n = P1_n*(1-prop) + prop*P3_n;
                end
            case 0
                P2_n = handles.all_handles{i}.new_position(2,:);
        end
        if handles.rectify_atual(i) == 1
            prop = handles.all_handles{i}.proportion;
            equi_coord_par = apply_transf_plot(handles.plot_pos{i},'Curved',handles.T_coefs(i,:),1);
            deriv = diff(equi_coord_par(:,1));
            index = find(abs(deriv) > 0.1);
            index = sort(index);
            if ~isempty(index)
                if (P1_n(1) - P3_n(1)) > P1_n(1)
                    P2_n = [P1_n(1)*(1-prop) + (prop)*(P3_n(1)+2*pi) P1_n(2)*(1-prop)+(prop)*P3_n(2)] ;
                elseif (P1_n(1) - P3_n(1))< P1_n(1)
                    P2_n = [P1_n(1)*(1-prop) + prop*(P3_n(1)-2*pi) P1_n(2)*(1-prop)+prop*P3_n(2)] ;
                end
            else
                P2_n = P1_n*(1-prop) + prop*P3_n;
            end
        end
        P1 = handles.all_handles{i}.position(1,:);
        P2 = handles.all_handles{i}.position(2,:);
        P3 = handles.all_handles{i}.position(3,:);
        handles.all_handles{i}.new_position = [P1_n;P2_n;P3_n];
        handles.T_coefs(i,:) = Calculate_Transformation([P1;P2;P3],[P1_n;P2_n;P3_n]);
        %%check if bone has a cluster
        if j ~= 3
            index_row = find(ismember(handles.point_index_clusters,[i j],'rows'));
            index_cluster = [];
            for k = 1:handles.num_of_clusters
                if find(handles.clustersInd{k} == index_row)
                    index_cluster = k;
                end
            end
            if ~isempty(index_cluster)
                for k = handles.clustersInd{index_cluster}
                    i_c_p = handles.point_index_clusters(k,:);
                    handles.all_handles{i_c_p(1)}.new_position(i_c_p(2),:) = handles.all_handles{i_c_p(1)}.position(i_c_p(2),:);
                    P1_n = handles.all_handles{i_c_p(1)}.new_position(1,:);
                    P3_n = handles.all_handles{i_c_p(1)}.new_position(3,:);
                    
                    switch get(handles.rectify,'Value')
                        case 1
                            prop = handles.all_handles{i_c_p(1)}.proportion;
                            equi_coord_par = apply_transf_plot(handles.plot_pos{i_c_p(1)},'Curved',handles.T_coefs(i_c_p(1),:),1);
                            deriv = diff(equi_coord_par(:,1));
                            index = find(abs(deriv) > 0.1);
                            index = sort(index);
                            if ~isempty(index)
                                if (P1_n(1) - P3_n(1))> P1_n(1)
                                    P2_n = [P1_n(1)*(1-prop) + prop*(P3_n(1)+2*pi) P1_n(2)*(1-prop)+prop*P3_n(2)] ;
                                elseif (P1_n(1) - P3_n(1))< P1_n(1)
                                    P2_n = [P1_n(1)*(1-prop) + prop*(P3_n(1)+2*pi) P1_n(2)*(1-prop)+prop*P3_n(2)] ;
                                end
                            else
                                P2_n = P1_n*(1-prop) + prop*P3_n;
                            end
                        case 0
                            if handles.rectify_atual(i_c_p(1)) == 1
                                prop = handles.all_handles{i_c_p(1)}.proportion;
                                equi_coord_par = apply_transf_plot(handles.plot_pos{i_c_p(1)},'Curved',handles.T_coefs(i_c_p(1),:),1);
                                deriv = diff(equi_coord_par(:,1));
                                index = find(abs(deriv) > 0.1);
                                index = sort(index);
                                if ~isempty(index)
                                    if (P1_n(1) - P3_n(1))> P1_n(1)
                                        P2_n = [P1_n(1)*(1-prop) + prop*(P3_n(1)+2*pi) P1_n(2)*(1-prop)+prop*P3_n(2)] ;
                                    elseif (P1_n(1) - P3_n(1))< P1_n(1)
                                        P2_n = [P1_n(1)*(1-prop) + prop*(P3_n(1)-2*pi) P1_n(2)*(1-prop)+prop*P3_n(2)] ;
                                    end
                                else
                                    P2_n = P1_n*(1-prop) + prop*P3_n;
                                end
                            else
                                P2_n = handles.all_handles{i_c_p(1)}.new_position(2,:);
                            end
                    end
                    
                    handles.all_handles{i_c_p(1)}.new_position = [P1_n; P2_n; P3_n];
                    handles.T_coefs(i_c_p(1),:) = Calculate_Transformation(handles.all_handles{i_c_p(1)}.position,handles.all_handles{i_c_p(1)}.new_position);
                end
            end
        end
    elseif isequal(handles.all_handles{i}.type,'Closed_Cage')
        P1_n = handles.all_handles{i}.position(1,:);
        P2_n = [P1_n(1)+0.3 P1_n(2)];
        P3_n = [P1_n(1)-0.3 P1_n(2)];
        P1 = handles.all_handles{i}.position(1,:);
        P2 = handles.all_handles{i}.position(2,:);
        P3 = handles.all_handles{i}.position(3,:);
        handles.all_handles{i}.new_position = [P1_n;P2_n;P3_n];
        handles.T_coefs(i,:) = Calculate_Transformation([P1;P2;P3],[P1_n;P2_n;P3_n]);
    end
    
    axes(handles.image)
    cla
    
    handles.matlabImage = biharmonic_moebius_sphere(handles.originalimage,handles.equi_points,handles.W,handles.T_coefs,...
        [handles.imagesize(1) handles.imagesize(2)]);
    
    handles.img_handle = image(handles.matlabImage);
    hold on
    handles.handles_plot = {};
    handles.cages_plot = {};
    for k = 1:length(handles.all_cages)
        for j = 1:length(handles.all_cages{k})-1
            index1 = handles.all_cages{k}(j);
            index2 = handles.all_cages{k}(j+1);
            Weight = handles.cage_weights{k}{j};
            transf = [handles.T_coefs(index1,:);handles.T_coefs(index2,:)];
            position = handles.cage_plot_pos{k}{j};
            equi_coord_par = apply_transf_plot(position,'Closed_Cage',transf,Weight);
            deriv = diff(equi_coord_par(:,1));
            index = find(abs(deriv) > 0.1);
            index = sort(index);
            equi_coord_par_pix = equi2pixels(equi_coord_par,handles.imagesize);
            if ~isempty(index)
                inicio = index(1);
                fim = index(length(index))+1;
                X1 = equi_coord_par_pix(1:inicio,1);
                Y1 = equi_coord_par_pix(1:inicio,2);
                X2 = equi_coord_par_pix(fim:100,1);
                Y2 = equi_coord_par_pix(fim:100,2);
                f = plot(X1,Y1,'Color',[1 1 1],'Parent', handles.image);
                f2 = plot(X2,Y2,'Color',[1 1 1],'Parent', handles.image);
                handles.cages_plot{end+1} = f;
                handles.cages_plot{end+1} = f2;
            else
                f = plot(equi_coord_par_pix(:,1), equi_coord_par_pix(:,2),'Color',[1 1 1],'Parent', handles.image);
                handles.cages_plot{end+1} = f;
            end
        end
    end
    for j= 1:handles.num_of_handles
        pix_pos = equi2pixels(handles.all_handles{j}.new_position,handles.imagesize);
        if isequal(handles.all_handles{j}.type,'Curved')
            equi_coord_par = apply_transf_plot(handles.plot_pos{j},'Curved',handles.T_coefs(j,:),1);
            equi_coord_par_pix = equi2pixels(equi_coord_par,handles.imagesize);
            deriv = diff(equi_coord_par(:,1));
            index = find(abs(deriv) > 0.1);
            index = sort(index);
            if ~isempty(index)
                inicio = index(1);
                fim = index(length(index))+1;
                X1 = equi_coord_par_pix(1:inicio,1);
                Y1 = equi_coord_par_pix(1:inicio,2);
                X2 = equi_coord_par_pix(fim:100,1);
                Y2 = equi_coord_par_pix(fim:100,2);
                l1 = plot(X1,Y1,'Color',[0.79 0 0.125],'Parent', handles.image);
                l2 = plot(X2,Y2,'Color',[0.79 0 0.125],'Parent', handles.image);
                f =  plot(pix_pos(1,1),pix_pos(1,2),'o','MarkerFaceColor',[0.79 0 0.125],'MarkerEdgeColor',[0.25 0.25 0.25],'MarkerSize',5,'Parent',handles.image);
                g =  plot(pix_pos(2,1),pix_pos(2,2),'o','MarkerFaceColor',[0.79 0 0.125],'MarkerEdgeColor',[0.25 0.25 0.25],'MarkerSize',5,'Parent',handles.image);
                h =  plot(pix_pos(3,1),pix_pos(3,2),'o','MarkerFaceColor',[0.79 0 0.125],'MarkerEdgeColor',[0.25 0.25 0.25],'MarkerSize',5,'Parent',handles.image);
                handles.handles_plot{end+1} = {f,g,h,l1,l2};
            else
                l = plot(equi_coord_par_pix(:,1),equi_coord_par_pix(:,2),'Color',[0.79 0 0.125]);
                f =  plot(pix_pos(1,1),pix_pos(1,2),'o','MarkerFaceColor',[0.79 0 0.125],'MarkerEdgeColor',[0.25 0.25 0.25],'MarkerSize',5,'Parent',handles.image);
                g =  plot(pix_pos(2,1),pix_pos(2,2),'o','MarkerFaceColor',[0.79 0 0.125],'MarkerEdgeColor',[0.25 0.25 0.25],'MarkerSize',5,'Parent',handles.image);
                h =  plot(pix_pos(3,1),pix_pos(3,2),'o','MarkerFaceColor',[0.79 0 0.125],'MarkerEdgeColor',[0.25 0.25 0.25],'MarkerSize',5,'Parent',handles.image);
                handles.handles_plot{end+1} = {f,g,h,l};
            end
            if handles.rectify_atual(j) == 1
                set(handles.handles_plot{j}{2},'Visible','off');
            end
            if ~handles.handles_visibility(j,1)
                set(handles.handles_plot{j}{1},'Visible','off');
            end
            if ~handles.handles_visibility(j,2)
                set(handles.handles_plot{j}{2},'Visible','off');
            end
            if ~handles.handles_visibility(j,3)
                set(handles.handles_plot{j}{3},'Visible','off');
            end
            set(handles.handles_plot{j}{1},'ButtonDownFcn',{@selectpoint,handles})
            set(handles.handles_plot{j}{2},'ButtonDownFcn',{@selectpoint,handles})
            set(handles.handles_plot{j}{3},'ButtonDownFcn',{@selectpoint,handles})
        elseif isequal(handles.all_handles{j}.type,'Point')
            f = arrow([pix_pos(1,1) pix_pos(1,2)],[pix_pos(2,1) pix_pos(2,2)],'Width',2,'Length',10,'FaceColor',[0.79 0 0.125],'EdgeColor',[0.25 0.25 0.25]);
            g = arrow([pix_pos(1,1) pix_pos(1,2)],[pix_pos(3,1) pix_pos(3,2)],'Width',2,'Length',10,'FaceColor',[0.79 0 0.125],'EdgeColor',[0.25 0.25 0.25]);
            h = plot(pix_pos(1,1),pix_pos(1,2),'o','MarkerFaceColor',[0.79 0 0.125],'MarkerEdgeColor',[0.25 0.25 0.25],'MarkerSize',5,'Parent',handles.image);
            handles.handles_plot{end+1} = {h,f,g};
            set(handles.handles_plot{j}{1},'UserData',[j 1]);
            set(handles.handles_plot{j}{2},'UserData',[j 2]);
            set(handles.handles_plot{j}{3},'UserData',[j 3]);
            set(handles.handles_plot{j}{1},'ButtonDownFcn',{@selectpoint,handles})
            set(handles.handles_plot{j}{2},'ButtonDownFcn',{@selectpoint,handles})
            set(handles.handles_plot{j}{3},'ButtonDownFcn',{@selectpoint,handles})
        elseif isequal(handles.all_handles{j}.type,'Closed_Cage')
            f = plot(pix_pos(1,1),pix_pos(1,2),'o','MarkerFaceColor',[0.79 0 0.125],'MarkerEdgeColor',[0.25 0.25 0.25],'MarkerSize',5,'Parent',handles.image);
            handles.handles_plot{j} = f;
            set(handles.handles_plot{j},'UserData',[j 1]);
            set(handles.handles_plot{j},'ButtonDownFcn',{@selectpoint,handles})
        end
    end
    axis off
    i = handles.current_handle(1);
    iii = handles.current_handle(2);
    if isequal(handles.all_handles{i}.type,'Point') && handles.current_handle(1,2) ~= 1
        set(handles.handles_plot{i}{iii},'FaceColor','y','EdgeColor',[0.25 0.25 0.25])
        set(handles.handles_plot{i}{iii}, 'ButtonDownFcn', {@image_ButtonDownFcn,handles});
    elseif isequal(handles.all_handles{i}.type,'Closed_Cage')
        set(handles.handles_plot{i},'MarkerFaceColor','y','MarkerEdgeColor',[0.25 0.25 0.25])
    else
        set(handles.handles_plot{i}{iii},'MarkerFaceColor','y','MarkerEdgeColor',[0.25 0.25 0.25])
    end
    hold off
    
end
setappdata(0,'Handles',handles);
hide_Callback(handles.hide,eventdata,handles);
Atualiza_handles();
setappdata(0,'Handles',handles);
guidata(hObject,handles);






% --- Executes on button press in save.
function save_Callback(hObject, eventdata, handles)
% hObject    handle to save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = getappdata(0,'Handles');
[filename, pathname] = uiputfile({'*.jpg'},'Save as');
image_to_save = biharmonic_moebius_sphere(handles.originalimage,handles.equi_points,handles.W,handles.T_coefs,...
    [handles.imagesize(1) handles.imagesize(2)]);
imwrite(image_to_save,strcat(pathname,filename),'JPEG');

% --- Executes on button press in rectify.
function rectify_Callback(hObject, eventdata, handles)
% hObject    handle to rectify (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = getappdata(0,'Handles');
% Hint: get(hObject,'Value') returns toggle state of rectify
switch get(hObject,'Value')
    case 0
        for i = 1:handles.num_of_handles
            if isequal(handles.all_handles{i}.type,'Curved') && handles.rectify_atual(i) == 0
                set(handles.handles_plot{i}{3},'Visible','on')
            end
        end
    case 1
        %%%weight for rectification
        w1 = 100;
        %%%weight for position
        w2 = 1;
        %%%weight for coupled endpoints
        w3 = 10;
        block_A = [];
        block_B = [];
        solution_index = {};
        solution = [];
        already_computed = [];
        for i = 1:handles.num_of_handles
            if isequal(handles.all_handles{i}.type,'Curved') && isempty(find(already_computed == i))
                %% bone without intersection
                P1_ = handles.all_handles{i}.new_position(1,:);
                P3_ = handles.all_handles{i}.new_position(3,:);
                t1_ = handles.all_handles{i}.proportion;
                dist = [sqrt((P1_(1)-P3_(1))^2) sqrt((P1_(1)-(P3_(1)+2*pi))^2) sqrt((P1_(1)-(P3_(1)-2*pi))^2)]; 
                [minimo,index] = min(dist);
                
                switch index
                    case 1
                    B = [w2*P1_(1); w2*P1_(2); 0; 0;w2*P3_(1); w2*P3_(2)];
                    case 2
                    B = [w2*P1_(1)-w1*2*pi*t1_*(1-t1_); w2*P1_(2); +w1*2*pi*t1_; 0;w2*P3_(1) - w1*2*pi*t1_*t1_; w2*P3_(2)];
                    case 3
                    B = [w2*P1_(1)+w1*2*pi*t1_*(1-t1_); w2*P1_(2); -w1*2*pi*t1_;  0   ;w2*P3_(1) + w1*2*pi*t1_*t1_; w2*P3_(2)];
                end
                A = [w1*(1-t1_)^2+w2       0       w1*(t1_-1)   0      w1*t1_*(1-t1_)      0       ; ...
                    0       w1*(1-t1_)^2+w2    0     w1*(t1_-1)      0       w1*t1_*(1-t1_); ...
                    w1*(t1_-1)         0           w1       0          -w1*t1_         0       ; ...
                    0          w1*(t1_-1)      0        w1           0          -w1*t1_    ; ...
                    w1*t1_*(1-t1_)        0        -w1*t1_     0       w1*t1_^2+w2        0       ; ...
                    0        w1*t1_*(1-t1_)    0     -w1*t1_         0        w1*t1_^2+w2] ;
                block_B = [block_B ; B];
                block_A = blkdiag(block_A,A);
                solution_index{end+1} = [i 1];
                solution_index{end+1} = [i 2];
                solution_index{end+1} = [i 3];
                already_computed = [already_computed i];
            elseif isequal(handles.all_handles{i}.type,'Point') && isempty(find(already_computed == i))
                %% point handle, only fixed position
                P1_ = handles.all_handles{i}.new_position(1,:);
                B = [w2*P1_(1); w2*P1_(2)];
                A = [w2 0 ; 0 w2];
                block_B = [block_B ; B];
                block_A = blkdiag(block_A,A);
                solution_index{end+1} = [i 1];
            elseif isequal(handles.all_handles{i}.type,'Closed_Cage') && isempty(find(already_computed == i))
                %% closed cage handle, only fixed position
                P1_ = handles.all_handles{i}.new_position(1,:);
                B = [w2*P1_(1); w2*P1_(2)];
                A = [w2 0 ; 0 w2];
                block_B = [block_B ; B];
                block_A = blkdiag(block_A,A);
                solution_index{end+1} = [i 1];
            end
        end
        for i = 1:handles.num_of_clusters
            %% coupled endpoints
            for j = 1:length(handles.clustersInd{i})-1
                for k = j+1:length(handles.clustersInd{i})
                    ind_j = handles.point_index_clusters(handles.clustersInd{i}(j),:);
                    ind_k = handles.point_index_clusters(handles.clustersInd{i}(k),:);
                    sol_ind_j = [];
                    sol_ind_k = [];
                    for l = 1:length(solution_index)
                        if ~isempty(intersect(solution_index{l},ind_j,'rows'))
                            sol_ind_j = l;
                            break
                        end
                    end
                    for l = 1:length(solution_index)
                        if ~isempty(intersect(solution_index{l},ind_k,'rows'))
                            sol_ind_k = l;
                            break
                        end
                    end
                    ind_i = 2*(sol_ind_j-1)+1;
                    ind_j = 2*(sol_ind_k-1)+1;
                    block_A(ind_i,ind_i) = block_A(ind_i,ind_i) + w3;
                    block_A(ind_i+1,ind_i+1) = block_A(ind_i+1,ind_i+1) + w3;
                    block_A(ind_i,ind_j) = block_A(ind_i,ind_j) - w3;
                    block_A(ind_i+1,ind_j+1) = block_A(ind_i+1,ind_j+1) - w3;
                    
                    block_A(ind_j,ind_j) = block_A(ind_j,ind_j) + w3;
                    block_A(ind_j+1,ind_j+1) = block_A(ind_j+1,ind_j+1) + w3;
                    block_A(ind_j,ind_i) = block_A(ind_j,ind_i) - w3;
                    block_A(ind_j+1,ind_i+1) = block_A(ind_j+1,ind_i+1) - w3;
                    
                end
            end
        end
        for i = 1:handles.num_of_handles
            for j = 1:handles.num_of_handles
                %% relative position of bones
                if isequal(handles.all_handles{i}.type,'Curved') && isequal(handles.all_handles{j}.type,'Curved') && i ~= j
                    P1_b = handles.all_handles{j}.new_position(1,:);
                    P3_b = handles.all_handles{j}.new_position(3,:);
                    bone_size = min([sqrt((P1_b(1) - P3_b(1))^2 + (P1_b(2) - P3_b(2))^2) ...
                                    sqrt((P1_b(1) - (P3_b(1)+2*pi))^2 + (P1_b(2) - P3_b(2))^2) ...
                                    sqrt((P1_b(1) - (P3_b(1)-2*pi))^2 + (P1_b(2) - P3_b(2))^2)]);
                    for k = 1:2
                        for l = 1:3
                           Pk_ = handles.all_handles{i}.new_position(2*(k-1)+1,:);
                           Pl_ = handles.all_handles{j}.new_position(l,:);
                           min_dist = [sqrt((Pk_(1) - Pl_(1))^2 + (Pk_(2) - Pl_(2))^2) ...
                                       sqrt((Pk_(1) - (Pl_(1)+2*pi))^2 + (Pk_(2) - Pl_(2))^2) ...
                                       sqrt((Pk_(1) - (Pl_(1)-2*pi))^2 + (Pk_(2) - Pl_(2))^2)];
                            minimo= min(min_dist);
                            if minimo < bone_size/2
                                ind_i = 0;
                                ind_j = 0;
                                
                                w4 = 5*kumaraswamy(minimo,1,5,bone_size/2,0);
                                for m = 1:length(solution_index)
                                    if ~isempty(intersect(solution_index{m},[i 2*(k-1)+1],'rows'))
                                        ind_i = m;
                                        break
                                    end
                                end
                                for m = 1:length(solution_index)
                                    if ~isempty(intersect(solution_index{m},[j l],'rows'))
                                        ind_j = m;
                                        break
                                    end
                                end
                                ind_i = 2*(ind_i-1)+1;
                                ind_j = 2*(ind_j-1)+1;
                                
                                block_A(ind_i,ind_i) = block_A(ind_i,ind_i) + w4;
                                block_A(ind_i+1,ind_i+1) = block_A(ind_i+1,ind_i+1) + w4;
                                block_A(ind_i,ind_j) = block_A(ind_i,ind_j) - w4;
                                block_A(ind_i+1,ind_j+1) = block_A(ind_i+1,ind_j+1) - w4;
                                
                                block_A(ind_j,ind_j) = block_A(ind_j,ind_j) + w4;
                                block_A(ind_j+1,ind_j+1) = block_A(ind_j+1,ind_j+1) + w4;
                                block_A(ind_j,ind_i) = block_A(ind_j,ind_i) - w4;
                                block_A(ind_j+1,ind_i+1) = block_A(ind_j+1,ind_i+1) - w4;
                                
                                block_B(ind_i) = block_B(ind_i) + w4*(handles.all_handles{i}.new_position(2*(k-1)+1,1) - handles.all_handles{j}.new_position(l,1));
                                block_B(ind_i+1) = block_B(ind_i+1) + w4*(handles.all_handles{i}.new_position(2*(k-1)+1,2) - handles.all_handles{j}.new_position(l,2));
                                block_B(ind_j) = block_B(ind_j) - w4*(handles.all_handles{i}.new_position(2*(k-1)+1,1) - handles.all_handles{j}.new_position(l,1));
                                block_B(ind_j+1) = block_B(ind_j+1) - w4*(handles.all_handles{i}.new_position(2*(k-1)+1,2) - handles.all_handles{j}.new_position(l,2));
                                                               
                            end
                        end
                    end
                elseif ( isequal(handles.all_handles{i}.type,'Point') || isequal(handles.all_handles{i}.type,'Closed_Cage')) && isequal(handles.all_handles{j}.type,'Curved') && i ~= j
                    %% relative position, point and cage handles
                    P1_b = handles.all_handles{j}.new_position(1,:);
                    P3_b = handles.all_handles{j}.new_position(3,:);
                    bone_size = min([sqrt((P1_b(1) - P3_b(1))^2 + (P1_b(2) - P3_b(2))^2) ...
                                    sqrt((P1_b(1) - (P3_b(1)+2*pi))^2 + (P1_b(2) - P3_b(2))^2) ...
                                    sqrt((P1_b(1) - (P3_b(1)-2*pi))^2 + (P1_b(2) - P3_b(2))^2)]);
                    for k = 1:3
                        Pk_ = handles.all_handles{i}.new_position(1,:);
                        Pl_ = handles.all_handles{j}.new_position(k,:);
                           min_dist = [sqrt((Pk_(1) - Pl_(1))^2 + (Pk_(2) - Pl_(2))^2) ...
                                       sqrt((Pk_(1) - (Pl_(1)+2*pi))^2 + (Pk_(2) - Pl_(2))^2) ...
                                       sqrt((Pk_(1) - (Pl_(1)-2*pi))^2 + (Pk_(2) - Pl_(2))^2)];
                           minimo= min(min_dist);
                        if minimo < bone_size/2
                            ind_i = 0;
                            ind_j = 0;

                            w4 = 5*kumaraswamy(minimo,1,5,bone_size/2,0);
                            for l = 1:length(solution_index)
                                if ~isempty(intersect(solution_index{l},[i 1],'rows'))
                                    ind_i = l;
                                    break
                                end
                            end
                            for l = 1:length(solution_index)
                                if ~isempty(intersect(solution_index{l},[j k],'rows'))
                                    ind_j = l;
                                    break
                                end
                            end
                            ind_i = 2*(ind_i-1)+1;
                            ind_j = 2*(ind_j-1)+1;
                            block_A(ind_i,ind_i) = block_A(ind_i,ind_i) + w4;
                            block_A(ind_i+1,ind_i+1) = block_A(ind_i+1,ind_i+1) + w4;
                            block_A(ind_i,ind_j) = block_A(ind_i,ind_j) - w4;
                            block_A(ind_i+1,ind_j+1) = block_A(ind_i+1,ind_j+1) - w4;
                            
                            block_A(ind_j,ind_j) = block_A(ind_j,ind_j) + w4;
                            block_A(ind_j+1,ind_j+1) = block_A(ind_j+1,ind_j+1) + w4;
                            block_A(ind_j,ind_i) = block_A(ind_j,ind_i) - w4;
                            block_A(ind_j+1,ind_i+1) = block_A(ind_j+1,ind_i+1) - w4;
                            
                            block_B(ind_i) = block_B(ind_i) + w4*(handles.all_handles{i}.new_position(1,1) - handles.all_handles{j}.new_position(k,1));
                            block_B(ind_i+1) = block_B(ind_i+1) + w4*(handles.all_handles{i}.new_position(1,2) - handles.all_handles{j}.new_position(k,2));
                            block_B(ind_j) = block_B(ind_j) - w4*(handles.all_handles{i}.new_position(1,1) - handles.all_handles{j}.new_position(k,1));
                            block_B(ind_j+1) = block_B(ind_j+1) - w4*(handles.all_handles{i}.new_position(1,2) - handles.all_handles{j}.new_position(k,2));
                        end
                    end
                end
            end
        end
        if ~isempty(block_B)
            solution = block_A\block_B;
        end
        if ~isempty(solution)
            solution = reshape(solution,[2,size(solution,1)*size(solution,2)/2])';
            for i = 1:length(solution_index)
                for j = 1:size(solution_index{i},1)
                    ind_1 = solution_index{i}(j,1);
                    ind_2 = solution_index{i}(j,2);
                    if (isequal(handles.all_handles{ind_1}.type,'Point') || isequal(handles.all_handles{ind_1}.type,'Closed_Cage'))
                        xlength = handles.all_handles{ind_1}.new_position(2,1) - handles.all_handles{ind_1}.new_position(1,1);
                        ylength = handles.all_handles{ind_1}.new_position(2,2) - handles.all_handles{ind_1}.new_position(1,2);
                        handles.all_handles{ind_1}.new_position(ind_2,:) = solution(i,:);
                        P1_n = handles.all_handles{ind_1}.new_position(1,:);
                        P2_n = [P1_n(1)+xlength P1_n(2)+ylength];
                        P3_n = [P1_n(1)-xlength P1_n(2)-ylength];
                        handles.all_handles{ind_1}.new_position(2,:) = P2_n;
                        handles.all_handles{ind_1}.new_position(3,:) = P3_n;
                    else
                        handles.all_handles{ind_1}.new_position(ind_2,:) = solution(i,:);
                    end
                end
            end
            for i = 1:handles.num_of_handles
                handles.T_coefs(i,:) = Calculate_Transformation(handles.all_handles{i}.position,handles.all_handles{i}.new_position);
            end
        end
        for i = 1:handles.num_of_handles
            if (isequal(handles.all_handles{i}.type,'Curved'))
                handles.all_handles{i}.new_mesh_position = apply_transf_plot(handles.all_handles{i}.mesh_position,'Curved',handles.T_coefs(i,:),1);
            end
        end
handles.matlabImage = biharmonic_moebius_sphere(handles.originalimage,handles.V,handles.F,handles.all_handles,handles.T_coefs,...
    [handles.imagesize(1) handles.imagesize(2)]);
        axes(handles.image)
        cla

        handles.img_handle = image(handles.matlabImage);
        hold on
        
        handles.handles_plot = {};
        handles.cages_plot = {};
        for k = 1:length(handles.all_cages)
            for j = 1:length(handles.all_cages{k})-1
                index1 = handles.all_cages{k}(j);
                index2 = handles.all_cages{k}(j+1);
                Weight = handles.cage_weights{k}{j};
                transf = [handles.T_coefs(index1,:);handles.T_coefs(index2,:)];
                position = handles.cage_plot_pos{k}{j};
                equi_coord_par = apply_transf_plot(position,'Closed_Cage',transf,Weight);
                deriv = diff(equi_coord_par(:,1));
                index = find(abs(deriv) > 0.1);
                index = sort(index);
                equi_coord_par_pix = equi2pixels(equi_coord_par,handles.imagesize);
                if ~isempty(index)
                    inicio = index(1);
                    fim = index(length(index))+1;
                    X1 = equi_coord_par_pix(1:inicio,1);
                    Y1 = equi_coord_par_pix(1:inicio,2);
                    X2 = equi_coord_par_pix(fim:100,1);
                    Y2 = equi_coord_par_pix(fim:100,2);
                    f = plot(X1,Y1,'Color',[1 1 1],'Parent', handles.image);
                    f2 = plot(X2,Y2,'Color',[1 1 1],'Parent', handles.image);
                    handles.cages_plot{end+1} = f;
                    handles.cages_plot{end+1} = f2;
                else
                    f = plot(equi_coord_par_pix(:,1), equi_coord_par_pix(:,2),'Color',[1 1 1],'Parent', handles.image);
                    handles.cages_plot{end+1} = f;
                end
            end
        end
        for j= 1:handles.num_of_handles
            pix_pos = equi2pixels(handles.all_handles{j}.new_position,handles.imagesize);
            if isequal(handles.all_handles{j}.type,'Curved')
                equi_coord_par = apply_transf_plot(handles.plot_pos{j},'Curved',handles.T_coefs(j,:),1);
                equi_coord_par_pix = equi2pixels(equi_coord_par,handles.imagesize);
                deriv = diff(equi_coord_par(:,1));
                index = find(abs(deriv) > 0.1);
                index = sort(index);
                if ~isempty(index)
                    inicio = index(1);
                    fim = index(length(index))+1;
                    X1 = equi_coord_par_pix(1:inicio,1);
                    Y1 = equi_coord_par_pix(1:inicio,2);
                    X2 = equi_coord_par_pix(fim:100,1);
                    Y2 = equi_coord_par_pix(fim:100,2);
                    l1 = plot(X1,Y1,'Color',[0.79 0 0.125],'Parent', handles.image);
                    l2 = plot(X2,Y2,'Color',[0.79 0 0.125],'Parent', handles.image);
                    f =  plot(pix_pos(1,1),pix_pos(1,2),'o','MarkerFaceColor',[0.79 0 0.125],'MarkerEdgeColor',[0.25 0.25 0.25],'MarkerSize',5,'Parent',handles.image);
                    g =  plot(pix_pos(2,1),pix_pos(2,2),'o','MarkerFaceColor',[0.79 0 0.125],'MarkerEdgeColor',[0.25 0.25 0.25],'MarkerSize',5,'Parent',handles.image);
                    h =  plot(pix_pos(3,1),pix_pos(3,2),'o','MarkerFaceColor',[0.79 0 0.125],'MarkerEdgeColor',[0.25 0.25 0.25],'MarkerSize',5,'Parent',handles.image);
                    handles.handles_plot{end+1} = {f,g,h,l1,l2};
                else
                    l = plot(equi_coord_par_pix(:,1),equi_coord_par_pix(:,2),'Color',[0.79 0 0.125]);
                    f =  plot(pix_pos(1,1),pix_pos(1,2),'o','MarkerFaceColor',[0.79 0 0.125],'MarkerEdgeColor',[0.25 0.25 0.25],'MarkerSize',5,'Parent',handles.image);
                    g =  plot(pix_pos(2,1),pix_pos(2,2),'o','MarkerFaceColor',[0.79 0 0.125],'MarkerEdgeColor',[0.25 0.25 0.25],'MarkerSize',5,'Parent',handles.image);
                    h =  plot(pix_pos(3,1),pix_pos(3,2),'o','MarkerFaceColor',[0.79 0 0.125],'MarkerEdgeColor',[0.25 0.25 0.25],'MarkerSize',5,'Parent',handles.image);
                    handles.handles_plot{end+1} = {f,g,h,l};
                end
                if handles.rectify_atual(j) == 1
                    set(handles.handles_plot{j}{2},'Visible','off');
                end
                if ~handles.handles_visibility(j,1)
                    set(handles.handles_plot{j}{1},'Visible','off');
                end
                if ~handles.handles_visibility(j,2)
                    set(handles.handles_plot{j}{2},'Visible','off');
                end
                if ~handles.handles_visibility(j,3)
                    set(handles.handles_plot{j}{3},'Visible','off');
                end
                set(handles.handles_plot{j}{1},'ButtonDownFcn',{@selectpoint,handles})
                set(handles.handles_plot{j}{2},'ButtonDownFcn',{@selectpoint,handles})
                set(handles.handles_plot{j}{3},'ButtonDownFcn',{@selectpoint,handles})
            elseif isequal(handles.all_handles{j}.type,'Point')
                f = arrow([pix_pos(1,1) pix_pos(1,2)],[pix_pos(2,1) pix_pos(2,2)],'Width',2,'Length',10,'FaceColor',[0.79 0 0.125],'EdgeColor',[0.25 0.25 0.25]);
                g = arrow([pix_pos(1,1) pix_pos(1,2)],[pix_pos(3,1) pix_pos(3,2)],'Width',2,'Length',10,'FaceColor',[0.79 0 0.125],'EdgeColor',[0.25 0.25 0.25]);
                h = plot(pix_pos(1,1),pix_pos(1,2),'o','MarkerFaceColor',[0.79 0 0.125],'MarkerEdgeColor',[0.25 0.25 0.25],'MarkerSize',5,'Parent',handles.image);
                handles.handles_plot{end+1} = {h,f,g};
                set(handles.handles_plot{j}{1},'UserData',[j 1]);
                set(handles.handles_plot{j}{2},'UserData',[j 2]);
                set(handles.handles_plot{j}{3},'UserData',[j 3]);
                set(handles.handles_plot{j}{1},'ButtonDownFcn',{@selectpoint,handles})
                set(handles.handles_plot{j}{2},'ButtonDownFcn',{@selectpoint,handles})
                set(handles.handles_plot{j}{3},'ButtonDownFcn',{@selectpoint,handles})
            elseif isequal(handles.all_handles{j}.type,'Closed_Cage')
                f = plot(pix_pos(1,1),pix_pos(1,2),'o','MarkerFaceColor',[0.79 0 0.125],'MarkerEdgeColor',[0.25 0.25 0.25],'MarkerSize',5,'Parent',handles.image);
                handles.handles_plot{j} = f;
                set(handles.handles_plot{j},'UserData',[j 1]);
                set(handles.handles_plot{j},'ButtonDownFcn',{@selectpoint,handles})
            end
        end
        axis off
        if ~isempty(handles.current_handle)
            i = handles.current_handle(1);
            iii = handles.current_handle(2);
            if isequal(handles.all_handles{i}.type,'Point') && handles.current_handle(1,2) ~= 1
                set(handles.handles_plot{i}{iii},'FaceColor','y','EdgeColor',[0.25 0.25 0.25])
                set(handles.handles_plot{i}{iii}, 'ButtonDownFcn', {@image_ButtonDownFcn,handles});
            elseif isequal(handles.all_handles{i}.type,'Closed_Cage')
                set(handles.handles_plot{i},'MarkerFaceColor','y','MarkerEdgeColor',[0.25 0.25 0.25])
            else
                set(handles.handles_plot{i}{iii},'MarkerFaceColor','y','MarkerEdgeColor',[0.25 0.25 0.25])
            end
        end
        hold off
end


setappdata(0,'Handles',handles);
hide_Callback(handles.hide,eventdata,handles);
Atualiza_handles();
setappdata(0,'Handles',handles);
guidata(hObject,handles);


% --- Executes on button press in hide.
function hide_Callback(hObject, eventdata, handles)
% hObject    handle to hide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
status = get(hObject,'Value');
handles = getappdata(0,'Handles');

switch status
    case 1
        for i = 1:handles.num_of_handles
            if ~isequal(handles.all_handles{i}.type,'Closed_Cage')
                for j = 1:size(handles.handles_plot{i},2)
                    set(handles.handles_plot{i}{j},'Visible','off')
                end
            else
                set(handles.handles_plot{i},'Visible','off')
            end
        end
        for i = 1:length(handles.cages_plot)
            set(handles.cages_plot{i},'Visible','off')
        end
    case 0
        for i = 1:handles.num_of_handles
            if ~isequal(handles.all_handles{i}.type,'Closed_Cage')
                for j = 1:size(handles.handles_plot{i},2)
                    set(handles.handles_plot{i}{j},'Visible','on')
                end
            else
                set(handles.handles_plot{i},'Visible','on')
            end
            if get(handles.rectify,'Value') == 1 && isequal(handles.all_handles{i}.type,'Curved')
                set(handles.handles_plot{i}{2},'Visible','off')
            elseif isequal(handles.all_handles{i}.type,'Curved') && handles.rectify_atual(i) == 1
                set(handles.handles_plot{i}{2},'Visible','off')
            end
        end
        for i = 1:length(handles.cages_plot)
            set(handles.cages_plot{i},'Visible','on')
        end
end


setappdata(0,'Handles',handles);
guidata(hObject,handles);
% Hint: get(hObject,'Value') returns toggle state of hide

% --- Outputs from this function are returned to the command line.
function varargout = transformation_interface_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function figure1_WindowButtonDownFcn(hObject,eventdata,handles)


% --- Executes on button press in rectify_current_line.
function rectify_current_line_Callback(hObject, eventdata, handles)
% hObject    handle to rectify_current_line (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
status = get(hObject,'Value');
handles = getappdata(0,'Handles');
if ~isempty(handles.current_handle)
    ind = handles.current_handle(1);
    ind2 = handles.current_handle(2);
    if isequal(handles.all_handles{ind}.type,'Curved')
        switch status
            case 1
                %%%weight for recfication
                w1 = 100;
                %%%weight for position
                w2 = 1;
                %%%weight for coupled endpoints
                w3 = 10;
                block_A = [];
                block_B = [];
                solution_index = {};
                solution = [];
                already_computed = [];
                handles.rectify_atual(ind) = 1;
                if ind2 == 2 || ind2 == 4
                    index_row = find(ismember(handles.point_index_clusters,[ind ind2-1],'rows'));
                    index_cluster = [];
                    for k = 1:handles.num_of_clusters
                        if find(handles.clustersInd{k} == index_row)
                            index_cluster = k;
                        end
                    end
                    if ~isempty(index_cluster)
                        for k = handles.clustersInd{index_cluster}
                            i_c_p = handles.point_index_clusters(k,:);
                            handles.rectify_atual(i_c_p(1)) = 1;
                        end
                    end
                end
                for i = 1:handles.num_of_handles
                    if isequal(handles.all_handles{i}.type,'Curved') && isempty(find(already_computed == i)) && handles.rectify_atual(i) == 1
                        %% bone rectification
                        P1_ = handles.all_handles{i}.new_position(1,:);
                        P3_ = handles.all_handles{i}.new_position(3,:);
                        t1_ = handles.all_handles{i}.proportion;
                        sphere_coord = equi2sphere(P1_);
                        sphere_coord2 = equi2sphere(P3_);
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
                        
                        
                        if isempty(index)
                            B = [w2*P1_(1); w2*P1_(2); 0; 0;w2*P3_(1); w2*P3_(2)];
                        elseif ~isempty(index) && (P1_(1)-P3_(1)) > P1_(1)
                            B = [w2*P1_(1)-w1*2*pi*t1_*(1-t1_); w2*P1_(2); +w1*2*pi*t1_; 0;w2*P3_(1) - w1*2*pi*t1_*t1_; w2*P3_(2)];
                        elseif ~isempty(index) && (P1_(1)-P3_(1)) < P1_(1)
                            B = [w2*P1_(1)+w1*2*pi*t1_*(1-t1_); w2*P1_(2); -w1*2*pi*t1_;  0   ;w2*P3_(1) + w1*2*pi*t1_*t1_; w2*P3_(2)];
                        end
                        A = [w1*(1-t1_)^2+w2       0       w1*(t1_-1)   0      w1*t1_*(1-t1_)      0       ; ...
                            0       w1*(1-t1_)^2+w2    0     w1*(t1_-1)      0       w1*t1_*(1-t1_); ...
                            w1*(t1_-1)         0           w1       0          -w1*t1_         0       ; ...
                            0          w1*(t1_-1)      0        w1           0          -w1*t1_    ; ...
                            w1*t1_*(1-t1_)        0        -w1*t1_     0       w1*t1_^2+w2        0       ; ...
                            0        w1*t1_*(1-t1_)    0     -w1*t1_         0        w1*t1_^2+w2] ;
                        block_B = [block_B ; B];
                        block_A = blkdiag(block_A,A);
                        solution_index{end+1} = [i 1];
                        solution_index{end+1} = [i 2];
                        solution_index{end+1} = [i 3];
                        already_computed = [already_computed i];
                    elseif isequal(handles.all_handles{i}.type,'Curved') && isempty(find(already_computed == i)) && handles.rectify_atual(i) == 0
                        %% non rectified bone, only fixed position
                        P1_ = handles.all_handles{i}.new_position(1,:);
                        P2_ = handles.all_handles{i}.new_position(2,:);
                        P3_ = handles.all_handles{i}.new_position(3,:);
                        
                        B = [w2*P1_(1); w2*P1_(2); w2*P2_(1); w2*P2_(2);w2*P3_(1); w2*P3_(2)];
                        A = diag([w2 w2 w2 w2 w2 w2]);
                        block_B = [block_B ; B];
                        block_A = blkdiag(block_A,A);
                        solution_index{end+1} = [i 1];
                        solution_index{end+1} = [i 2];
                        solution_index{end+1} = [i 3];
                        already_computed = [already_computed i];
                    elseif isequal(handles.all_handles{i}.type,'Point') && isempty(find(already_computed == i))
                        %% point handle only fixed position
                        P1_ = handles.all_handles{i}.new_position(1,:);
                        B = [w2*P1_(1); w2*P1_(2)];
                        A = [w2 0 ; 0 w2];
                        block_B = [block_B ; B];
                        block_A = blkdiag(block_A,A);
                        solution_index{end+1} = [i 1];
                    elseif isequal(handles.all_handles{i}.type,'Closed_Cage') && isempty(find(already_computed == i))
                        %% closed cage handle only fixed position
                        P1_ = handles.all_handles{i}.new_position(1,:);
                        B = [w2*P1_(1); w2*P1_(2)];
                        A = [w2 0 ; 0 w2];
                        block_B = [block_B ; B];
                        block_A = blkdiag(block_A,A);
                        solution_index{end+1} = [i 1];
                    end
                end
                for i = 1:handles.num_of_clusters
                    %% coupled endpoints
                    for j = 1:length(handles.clustersInd{i})-1
                        for k = j+1:length(handles.clustersInd{i})
                            ind_j = handles.point_index_clusters(handles.clustersInd{i}(j),:);
                            ind_k = handles.point_index_clusters(handles.clustersInd{i}(k),:);
                            sol_ind_j = [];
                            sol_ind_k = [];
                            for l = 1:length(solution_index)
                                if ~isempty(intersect(solution_index{l},ind_j,'rows'))
                                    sol_ind_j = l;
                                    break
                                end
                            end
                            for l = 1:length(solution_index)
                                if ~isempty(intersect(solution_index{l},ind_k,'rows'))
                                    sol_ind_k = l;
                                    break
                                end
                            end
                            ind_i = 2*(sol_ind_j-1)+1;
                            ind_j = 2*(sol_ind_k-1)+1;
                            block_A(ind_i,ind_i) = block_A(ind_i,ind_i) + w3;
                            block_A(ind_i+1,ind_i+1) = block_A(ind_i+1,ind_i+1) + w3;
                            block_A(ind_i,ind_j) = block_A(ind_i,ind_j) - w3;
                            block_A(ind_i+1,ind_j+1) = block_A(ind_i+1,ind_j+1) - w3;
                            
                            block_A(ind_j,ind_j) = block_A(ind_j,ind_j) + w3;
                            block_A(ind_j+1,ind_j+1) = block_A(ind_j+1,ind_j+1) + w3;
                            block_A(ind_j,ind_i) = block_A(ind_j,ind_i) - w3;
                            block_A(ind_j+1,ind_i+1) = block_A(ind_j+1,ind_i+1) - w3;
                            
                        end
                    end
                end
                for i = 1:handles.num_of_handles
                    for j = 1:handles.num_of_handles
                        %% relative position
                        if isequal(handles.all_handles{i}.type,'Curved') && isequal(handles.all_handles{j}.type,'Curved') && i ~= j && handles.rectify_atual(i) == 1
                            dist_11 = sqrt(sum((handles.all_handles{i}.new_position(1,:) - handles.all_handles{j}.new_position(1,:)).^2));
                            dist_12 = sqrt(sum((handles.all_handles{i}.new_position(1,:) - handles.all_handles{j}.new_position(2,:)).^2));
                            dist_13 = sqrt(sum((handles.all_handles{i}.new_position(1,:) - handles.all_handles{j}.new_position(3,:)).^2));
                            dist_31 = sqrt(sum((handles.all_handles{i}.new_position(3,:) - handles.all_handles{j}.new_position(1,:)).^2));
                            dist_32 = sqrt(sum((handles.all_handles{i}.new_position(3,:) - handles.all_handles{j}.new_position(2,:)).^2));
                            dist_33 = sqrt(sum((handles.all_handles{i}.new_position(3,:) - handles.all_handles{j}.new_position(3,:)).^2));
                            bone_size = sqrt(sum((handles.all_handles{j}.new_position(1,:) - handles.all_handles{j}.new_position(3,:)).^2));
                            dist = [ dist_11 dist_12 dist_13; dist_31 dist_32 dist_33];
                            for k = 1:2
                                for l = 1:3
                                    if dist(k,l) < bone_size/2
                                        ind_i = 0;
                                        ind_j = 0;

                                        w4 = 5*kumaraswamy(dist(k,l),1,5,bone_size/2,0);
                                        for m = 1:length(solution_index)
                                            if ~isempty(intersect(solution_index{m},[i 2*(k-1)+1],'rows'))
                                                ind_i = m;
                                                break
                                            end
                                        end
                                        for m = 1:length(solution_index)
                                            if ~isempty(intersect(solution_index{m},[j l],'rows'))
                                                ind_j = m;
                                                break
                                            end
                                        end
                                        ind_i = 2*(ind_i-1)+1;
                                        ind_j = 2*(ind_j-1)+1;
                                        
                                        block_A(ind_i,ind_i) = block_A(ind_i,ind_i) + w4;
                                        block_A(ind_i+1,ind_i+1) = block_A(ind_i+1,ind_i+1) + w4;
                                        block_A(ind_i,ind_j) = block_A(ind_i,ind_j) - w4;
                                        block_A(ind_i+1,ind_j+1) = block_A(ind_i+1,ind_j+1) - w4;
                                        
                                        block_A(ind_j,ind_j) = block_A(ind_j,ind_j) + w4;
                                        block_A(ind_j+1,ind_j+1) = block_A(ind_j+1,ind_j+1) + w4;
                                        block_A(ind_j,ind_i) = block_A(ind_j,ind_i) - w4;
                                        block_A(ind_j+1,ind_i+1) = block_A(ind_j+1,ind_i+1) - w4;
                                        
                                        block_B(ind_i) = block_B(ind_i) + w4*(handles.all_handles{i}.new_position(2*(k-1)+1,1) - handles.all_handles{j}.new_position(l,1));
                                        block_B(ind_i+1) = block_B(ind_i+1) + w4*(handles.all_handles{i}.new_position(2*(k-1)+1,2) - handles.all_handles{j}.new_position(l,2));
                                        block_B(ind_j) = block_B(ind_j) - w4*(handles.all_handles{i}.new_position(2*(k-1)+1,1) - handles.all_handles{j}.new_position(l,1));
                                        block_B(ind_j+1) = block_B(ind_j+1) - w4*(handles.all_handles{i}.new_position(2*(k-1)+1,2) - handles.all_handles{j}.new_position(l,2));
                                    end
                                end
                            end
                            
                        elseif isequal(handles.all_handles{i}.type,'Curved') && isequal(handles.all_handles{j}.type,'Curved') && i ~= j && handles.rectify_atual(i) == 0
                            %% relative position
                            dist_11 = sqrt(sum((handles.all_handles{i}.new_position(1,:) - handles.all_handles{j}.new_position(1,:)).^2));
                            dist_12 = sqrt(sum((handles.all_handles{i}.new_position(1,:) - handles.all_handles{j}.new_position(2,:)).^2));
                            dist_13 = sqrt(sum((handles.all_handles{i}.new_position(1,:) - handles.all_handles{j}.new_position(3,:)).^2));
                            dist_21 = sqrt(sum((handles.all_handles{i}.new_position(2,:) - handles.all_handles{j}.new_position(1,:)).^2));
                            dist_22 = sqrt(sum((handles.all_handles{i}.new_position(2,:) - handles.all_handles{j}.new_position(2,:)).^2));
                            dist_23 = sqrt(sum((handles.all_handles{i}.new_position(2,:) - handles.all_handles{j}.new_position(3,:)).^2));
                            dist_31 = sqrt(sum((handles.all_handles{i}.new_position(3,:) - handles.all_handles{j}.new_position(1,:)).^2));
                            dist_32 = sqrt(sum((handles.all_handles{i}.new_position(3,:) - handles.all_handles{j}.new_position(2,:)).^2));
                            dist_33 = sqrt(sum((handles.all_handles{i}.new_position(3,:) - handles.all_handles{j}.new_position(3,:)).^2));
                            bone_size = sqrt(sum((handles.all_handles{j}.new_position(1,:) - handles.all_handles{j}.new_position(3,:)).^2));
                            dist = [ dist_11 dist_12 dist_13; dist_21 dist_22 dist_23; dist_31 dist_32 dist_33];
                            for k = 1:3
                                for l = 1:3
                                    if dist(k,l) < bone_size/2
                                        ind_i = 0;
                                        ind_j = 0;

                                        w4 = 5*kumaraswamy(dist(k,l),1,5,bone_size/2,0);
                                        for m = 1:length(solution_index)
                                            if ~isempty(intersect(solution_index{m},[i k],'rows'))
                                                ind_i = m;
                                                break
                                            end
                                        end
                                        for m = 1:length(solution_index)
                                            if ~isempty(intersect(solution_index{m},[j l],'rows'))
                                                ind_j = m;
                                                break
                                            end
                                        end
                                        ind_i = 2*(ind_i-1)+1;
                                        ind_j = 2*(ind_j-1)+1;
                                        
                                        block_A(ind_i,ind_i) = block_A(ind_i,ind_i) + w4;
                                        block_A(ind_i+1,ind_i+1) = block_A(ind_i+1,ind_i+1) + w4;
                                        block_A(ind_i,ind_j) = block_A(ind_i,ind_j) - w4;
                                        block_A(ind_i+1,ind_j+1) = block_A(ind_i+1,ind_j+1) - w4;
                                        
                                        block_A(ind_j,ind_j) = block_A(ind_j,ind_j) + w4;
                                        block_A(ind_j+1,ind_j+1) = block_A(ind_j+1,ind_j+1) + w4;
                                        block_A(ind_j,ind_i) = block_A(ind_j,ind_i) - w4;
                                        block_A(ind_j+1,ind_i+1) = block_A(ind_j+1,ind_i+1) - w4;
                                        
                                        block_B(ind_i) = block_B(ind_i) + w4*(handles.all_handles{i}.new_position(k,1) - handles.all_handles{j}.new_position(l,1));
                                        block_B(ind_i+1) = block_B(ind_i+1) + w4*(handles.all_handles{i}.new_position(k,2) - handles.all_handles{j}.new_position(l,2));
                                        block_B(ind_j) = block_B(ind_j) - w4*(handles.all_handles{i}.new_position(k,1) - handles.all_handles{j}.new_position(l,1));
                                        block_B(ind_j+1) = block_B(ind_j+1) - w4*(handles.all_handles{i}.new_position(k,2) - handles.all_handles{j}.new_position(l,2));
                                    end
                                end
                            end
                        elseif ( isequal(handles.all_handles{i}.type,'Point') || isequal(handles.all_handles{i}.type,'Closed_Cage')) && isequal(handles.all_handles{j}.type,'Curved') && i ~= j
                            %% relative position
                            dist_11 = sqrt(sum((handles.all_handles{i}.new_position(1,:) - handles.all_handles{j}.new_position(1,:)).^2));
                            dist_12 = sqrt(sum((handles.all_handles{i}.new_position(1,:) - handles.all_handles{j}.new_position(2,:)).^2));
                            dist_13 = sqrt(sum((handles.all_handles{i}.new_position(1,:) - handles.all_handles{j}.new_position(3,:)).^2));
                            bone_size = sqrt(sum((handles.all_handles{j}.new_position(1,:) - handles.all_handles{j}.new_position(3,:)).^2));
                            dist = [ dist_11 dist_12 dist_13];
                            for k = 1:3
                                if dist(k) < bone_size/2
                                    ind_i = 0;
                                    ind_j = 0;

                                    w4 = 5*kumaraswamy(dist(k),1,5,bone_size/2,0);
                                    for l = 1:length(solution_index)
                                        if ~isempty(intersect(solution_index{l},[i 1],'rows'))
                                            ind_i = l;
                                            break
                                        end
                                    end
                                    for l = 1:length(solution_index)
                                        if ~isempty(intersect(solution_index{l},[j k],'rows'))
                                            ind_j = l;
                                            break
                                        end
                                    end
                                    ind_i = 2*(ind_i-1)+1;
                                    ind_j = 2*(ind_j-1)+1;
                                    block_A(ind_i,ind_i) = block_A(ind_i,ind_i) + w4;
                                    block_A(ind_i+1,ind_i+1) = block_A(ind_i+1,ind_i+1) + w4;
                                    block_A(ind_i,ind_j) = block_A(ind_i,ind_j) - w4;
                                    block_A(ind_i+1,ind_j+1) = block_A(ind_i+1,ind_j+1) - w4;
                                    
                                    block_A(ind_j,ind_j) = block_A(ind_j,ind_j) + w4;
                                    block_A(ind_j+1,ind_j+1) = block_A(ind_j+1,ind_j+1) + w4;
                                    block_A(ind_j,ind_i) = block_A(ind_j,ind_i) - w4;
                                    block_A(ind_j+1,ind_i+1) = block_A(ind_j+1,ind_i+1) - w4;
                                    
                                    block_B(ind_i) = block_B(ind_i) + w4*(handles.all_handles{i}.new_position(1,1) - handles.all_handles{j}.new_position(k,1));
                                    block_B(ind_i+1) = block_B(ind_i+1) + w4*(handles.all_handles{i}.new_position(1,2) - handles.all_handles{j}.new_position(k,2));
                                    block_B(ind_j) = block_B(ind_j) - w4*(handles.all_handles{i}.new_position(1,1) - handles.all_handles{j}.new_position(k,1));
                                    block_B(ind_j+1) = block_B(ind_j+1) - w4*(handles.all_handles{i}.new_position(1,2) - handles.all_handles{j}.new_position(k,2));
                                end
                            end
                        end
                    end
                end
                if ~isempty(block_B)
                    solution = block_A\block_B;
                end
                if ~isempty(solution)
                    solution = reshape(solution,[2,size(solution,1)*size(solution,2)/2])';
                    for i = 1:length(solution_index)
                        for j = 1:size(solution_index{i},1)
                            ind_1 = solution_index{i}(j,1);
                            ind_2 = solution_index{i}(j,2);
                            if (isequal(handles.all_handles{ind_1}.type,'Point') || isequal(handles.all_handles{ind_1}.type,'Closed_Cage'))
                                xlength = handles.all_handles{ind_1}.new_position(2,1) - handles.all_handles{ind_1}.new_position(1,1);
                                ylength = handles.all_handles{ind_1}.new_position(2,2) - handles.all_handles{ind_1}.new_position(1,2);
                                handles.all_handles{ind_1}.new_position(ind_2,:) = solution(i,:);
                                P1_n = handles.all_handles{ind_1}.new_position(1,:);
                                P2_n = [P1_n(1)+xlength P1_n(2)+ylength];
                                P3_n = [P1_n(1)-xlength P1_n(2)-ylength];
                                handles.all_handles{ind_1}.new_position(2,:) = P2_n;
                                handles.all_handles{ind_1}.new_position(3,:) = P3_n;
                            else
                                handles.all_handles{ind_1}.new_position(ind_2,:) = solution(i,:);
                            end
                        end
                    end
                    for i = 1:handles.num_of_handles
                        for j = 1:3
                            if handles.all_handles{i}.new_position(j,1) > pi
                               handles.all_handles{i}.new_position(j,1) = handles.all_handles{i}.new_position(j,1) - 2*pi;
                            elseif handles.all_handles{i}.new_position(j,1) < -pi
                                handles.all_handles{i}.new_position(j,1) = handles.all_handles{i}.new_position(j,1) + 2*pi;
                            end
                        end
                        handles.T_coefs(i,:) = Calculate_Transformation(handles.all_handles{i}.position,handles.all_handles{i}.new_position);
                    end
                end
                for i = 1:handles.num_of_handles
                    if (isequal(handles.all_handles{i}.type,'Curved'))
                        handles.all_handles{i}.new_mesh_position = apply_transf_plot(handles.all_handles{i}.mesh_position,'Curved',handles.T_coefs(i,:),1);
                    end
                end
                handles.matlabImage = biharmonic_moebius_sphere(handles.originalimage,handles.V,handles.F,handles.all_handles,handles.T_coefs,...
                [handles.imagesize(1) handles.imagesize(2)]);
                axes(handles.image)
                cla
                handles.img_handle = image(handles.matlabImage);
                hold on
                
                handles.handles_plot = {};
                handles.cages_plot = {};
                for k = 1:length(handles.all_cages)
                    for j = 1:length(handles.all_cages{k})-1
                        index1 = handles.all_cages{k}(j);
                        index2 = handles.all_cages{k}(j+1);
                        Weight = handles.cage_weights{k}{j};
                        transf = [handles.T_coefs(index1,:);handles.T_coefs(index2,:)];
                        position = handles.cage_plot_pos{k}{j};
                        equi_coord_par = apply_transf_plot(position,'Closed_Cage',transf,Weight);
                        deriv = diff(equi_coord_par(:,1));
                        index = find(abs(deriv) > 0.1);
                        index = sort(index);
                        equi_coord_par_pix = equi2pixels(equi_coord_par,handles.imagesize);
                        if ~isempty(index)
                            inicio = index(1);
                            fim = index(length(index))+1;
                            X1 = equi_coord_par_pix(1:inicio,1);
                            Y1 = equi_coord_par_pix(1:inicio,2);
                            X2 = equi_coord_par_pix(fim:100,1);
                            Y2 = equi_coord_par_pix(fim:100,2);
                            f = plot(X1,Y1,'Color',[1 1 1],'Parent', handles.image);
                            f2 = plot(X2,Y2,'Color',[1 1 1],'Parent', handles.image);
                            handles.cages_plot{end+1} = f;
                            handles.cages_plot{end+1} = f2;
                        else
                            f = plot(equi_coord_par_pix(:,1), equi_coord_par_pix(:,2),'Color',[1 1 1],'Parent', handles.image);
                            handles.cages_plot{end+1} = f;
                        end
                    end
                end
                for j= 1:handles.num_of_handles
                    pix_pos = equi2pixels(handles.all_handles{j}.new_position,handles.imagesize);
                    if isequal(handles.all_handles{j}.type,'Curved')
                        equi_coord_par = apply_transf_plot(handles.plot_pos{j},'Curved',handles.T_coefs(j,:),1);
                        equi_coord_par_pix = equi2pixels(equi_coord_par,handles.imagesize);
                        deriv = diff(equi_coord_par(:,1));
                        index = find(abs(deriv) > 0.1);
                        index = sort(index);
                        if ~isempty(index)
                            inicio = index(1);
                            fim = index(length(index))+1;
                            X1 = equi_coord_par_pix(1:inicio,1);
                            Y1 = equi_coord_par_pix(1:inicio,2);
                            X2 = equi_coord_par_pix(fim:100,1);
                            Y2 = equi_coord_par_pix(fim:100,2);
                            l1 = plot(X1,Y1,'Color',[0.79 0 0.125],'Parent', handles.image);
                            l2 = plot(X2,Y2,'Color',[0.79 0 0.125],'Parent', handles.image);
                            f =  plot(pix_pos(1,1),pix_pos(1,2),'o','MarkerFaceColor',[0.79 0 0.125],'MarkerEdgeColor',[0.25 0.25 0.25],'MarkerSize',5,'Parent',handles.image);
                            g =  plot(pix_pos(2,1),pix_pos(2,2),'o','MarkerFaceColor',[0.79 0 0.125],'MarkerEdgeColor',[0.25 0.25 0.25],'MarkerSize',5,'Parent',handles.image);
                            h =  plot(pix_pos(3,1),pix_pos(3,2),'o','MarkerFaceColor',[0.79 0 0.125],'MarkerEdgeColor',[0.25 0.25 0.25],'MarkerSize',5,'Parent',handles.image);
                            handles.handles_plot{end+1} = {f,g,h,l1,l2};
                        else
                            l = plot(equi_coord_par_pix(:,1),equi_coord_par_pix(:,2),'Color',[0.79 0 0.125]);
                            f =  plot(pix_pos(1,1),pix_pos(1,2),'o','MarkerFaceColor',[0.79 0 0.125],'MarkerEdgeColor',[0.25 0.25 0.25],'MarkerSize',5,'Parent',handles.image);
                            g =  plot(pix_pos(2,1),pix_pos(2,2),'o','MarkerFaceColor',[0.79 0 0.125],'MarkerEdgeColor',[0.25 0.25 0.25],'MarkerSize',5,'Parent',handles.image);
                            h =  plot(pix_pos(3,1),pix_pos(3,2),'o','MarkerFaceColor',[0.79 0 0.125],'MarkerEdgeColor',[0.25 0.25 0.25],'MarkerSize',5,'Parent',handles.image);
                            handles.handles_plot{end+1} = {f,g,h,l};
                        end
                        if handles.rectify_atual(j) == 1
                            set(handles.handles_plot{j}{2},'Visible','off');
                        end
                        if ~handles.handles_visibility(j,1)
                            set(handles.handles_plot{j}{1},'Visible','off');
                        end
                        if ~handles.handles_visibility(j,2)
                            set(handles.handles_plot{j}{2},'Visible','off');
                        end
                        if ~handles.handles_visibility(j,3)
                            set(handles.handles_plot{j}{3},'Visible','off');
                        end
                        set(handles.handles_plot{j}{1},'ButtonDownFcn',{@selectpoint,handles})
                        set(handles.handles_plot{j}{2},'ButtonDownFcn',{@selectpoint,handles})
                        set(handles.handles_plot{j}{3},'ButtonDownFcn',{@selectpoint,handles})
                    elseif isequal(handles.all_handles{j}.type,'Point')
                        f = arrow([pix_pos(1,1) pix_pos(1,2)],[pix_pos(2,1) pix_pos(2,2)],'Width',2,'Length',10,'FaceColor',[0.79 0 0.125],'EdgeColor',[0.25 0.25 0.25]);
                        g = arrow([pix_pos(1,1) pix_pos(1,2)],[pix_pos(3,1) pix_pos(3,2)],'Width',2,'Length',10,'FaceColor',[0.79 0 0.125],'EdgeColor',[0.25 0.25 0.25]);
                        h = plot(pix_pos(1,1),pix_pos(1,2),'o','MarkerFaceColor',[0.79 0 0.125],'MarkerEdgeColor',[0.25 0.25 0.25],'MarkerSize',5,'Parent',handles.image);
                        handles.handles_plot{end+1} = {h,f,g};
                        set(handles.handles_plot{j}{1},'UserData',[j 1]);
                        set(handles.handles_plot{j}{2},'UserData',[j 2]);
                        set(handles.handles_plot{j}{3},'UserData',[j 3]);
                        set(handles.handles_plot{j}{1},'ButtonDownFcn',{@selectpoint,handles})
                        set(handles.handles_plot{j}{2},'ButtonDownFcn',{@selectpoint,handles})
                        set(handles.handles_plot{j}{3},'ButtonDownFcn',{@selectpoint,handles})
                    elseif isequal(handles.all_handles{j}.type,'Closed_Cage')
                        f = plot(pix_pos(1,1),pix_pos(1,2),'o','MarkerFaceColor',[0.79 0 0.125],'MarkerEdgeColor',[0.25 0.25 0.25],'MarkerSize',5,'Parent',handles.image);
                        handles.handles_plot{j} = f;
                        set(handles.handles_plot{j},'UserData',[j 1]);
                        set(handles.handles_plot{j},'ButtonDownFcn',{@selectpoint,handles})
                    end
                end
                axis off
                i = handles.current_handle(1);
                iii = handles.current_handle(2);
                if isequal(handles.all_handles{i}.type,'Point') && handles.current_handle(1,2) ~= 1
                    set(handles.handles_plot{i}{iii},'FaceColor','y','EdgeColor',[0.25 0.25 0.25])
                    set(handles.handles_plot{i}{iii}, 'ButtonDownFcn', {@image_ButtonDownFcn,handles});
                elseif isequal(handles.all_handles{i}.type,'Closed_Cage')
                    set(handles.handles_plot{i},'MarkerFaceColor','y','MarkerEdgeColor',[0.25 0.25 0.25])
                else
                    set(handles.handles_plot{i}{iii},'MarkerFaceColor','y','MarkerEdgeColor',[0.25 0.25 0.25])
                end
                hold off
                
                
            case 0
                set(handles.handles_plot{ind}{2},'Visible','on')
                handles.rectify_atual(ind) = 0;
                if ind2 == 1 || ind2 == 3
                    index_row = find(ismember(handles.point_index_clusters,[ind ind2],'rows'));
                    index_cluster = [];
                    for k = 1:handles.num_of_clusters
                        if find(handles.clustersInd{k} == index_row)
                            index_cluster = k;
                        end
                    end
                    if ~isempty(index_cluster)
                        for k = handles.clustersInd{index_cluster}
                            i_c_p = handles.point_index_clusters(k,:);
                            handles.rectify_atual(i_c_p(1)) = 0;
                            set(handles.handles_plot{i_c_p(1)}{2},'Visible','on')
                        end
                    end
                end
        end
    else
        set(hObject,'Value',0);
    end
else
    set(hObject,'Value',0);
end
setappdata(0,'Handles',handles);
Atualiza_handles();
setappdata(0,'Handles',handles);
guidata(hObject,handles);

% Hint: get(hObject,'Value') returns toggle state of rectify_current_line

% --- Executes when figure1 is resized.
function figure1_SizeChangedFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

window_values = get(handles.figure1,'Position');

%  image
Pos_image = get(handles.image,'Position');
Pos_image(3) = min([2*(window_values(4)-80) window_values(3)-40]);
Pos_image(4) = (Pos_image(3)/2);
Pos_image(2) = round((window_values(4)-60)/2 - Pos_image(4)/2);
Pos_image(1) = (window_values(3)-Pos_image(3))/2;
set(handles.image,'Position',Pos_image);

Pos_undo = get(handles.undo,'Position');
Pos_undo(1) = 40;
Pos_undo(2) = window_values(4) - 60;
set(handles.undo,'Position',Pos_undo);

Pos_rectify = get(handles.rectify,'Position');
Pos_rectify(1) = 200;
Pos_rectify(2) = window_values(4) - 60;
set(handles.rectify,'Position',Pos_rectify);

Pos_hide = get(handles.hide,'Position');
Pos_hide(1) = 320;
Pos_hide(2) = window_values(4) - 60;
set(handles.hide,'Position',Pos_hide);

Pos_save = get(handles.save,'Position');
Pos_save(1) = window_values(3) - 160;
Pos_save(2) = window_values(4) - 60;
set(handles.save,'Position',Pos_save);

pos_rectify_atual = get(handles.rectify_current_line,'Position');
pos_rectify_atual(1) = 460;
pos_rectify_atual(2) = window_values(4) - 60;
set(handles.rectify_current_line,'Position',pos_rectify_atual);
