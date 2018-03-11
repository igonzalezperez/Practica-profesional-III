
function varargout = gui(varargin)
% GUI MATLAB code for GUI.fig
%      GUI, by itself, creates a new GUI or raises the existing
%      singleton*.
%
%      H = GUI returns the handle to a new GUI or the handle to
%      the existing singleton*.
%
%      GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI.M with the given input arguments.
%
%      GUI('Property','Value',...) creates a new GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI

% Last Modified by GUIDE v2.5 02-Feb-2018 12:44:49

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_OutputFcn, ...
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
function GUI_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI (see VARARGIN)

% Choose default command line output for GUI

mainfolder=cd;
export='Overallyangulos';
handles.export=export;
handles.mainfolder=mainfolder;
datafolder=[mainfolder '\Datasets'];
handles.datafolder=datafolder;
cd(mainfolder)
cd (datafolder)
handles.output = hObject;
set(handles.directorio_principal,'String',mainfolder)
set(handles.directorio_datos,'String',datafolder)
set(handles.dataexport,'String',export)
axes(handles.axes1)
delete( setdiff( findall(0, 'type', 'figure'), hObject ) );
guidata(hObject, handles);
% if strcmp(handles.datafolder,get(handles.directorio_datos,'String'))
%     cd(get(handles.directorio_datos,'String'))
% end
if nargin == 3
    initial_dir = pwd;
elseif nargin > 4
    if strcmpi(varargin{1},'dir')
        if exist(varargin{2},'dir')
            initial_dir = varargin{2};
        else
            errordlg('Input argument must be a valid directory','Input Argument Error!')
            return
        end
    else
        errordlg('Unrecognized input argument','Input Argument Error!');
        return;
    end
end
load_listbox(initial_dir,handles)




function varargout = GUI_OutputFcn(~, ~, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



%%%LISTBOX%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on double click press in lista_datos.
function lista_datos_Callback(~, ~, handles)
% hObject    handle to lista_datos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns lista_datos contents as cell array
%        contents{get(hObject,'Value')} returns selected item from lista_datos
mainfolder=handles.mainfolder;
datafolder=handles.datafolder;
get(handles.figure1,'SelectionType');
cd (datafolder)
if strcmp(get(handles.figure1,'SelectionType'),'open')
    index_selected = get(handles.lista_datos,'Value');
    file_list = get(handles.lista_datos,'String');
    filename = file_list{index_selected};
    handles.filename=filename;
    
    if length(filename)<4
        filename=[filename '...'];
    end
    
    if ~strcmp(filename(end-3:end),'.csv')
        warndlg('Seleccione un archivo compatible (.csv)');
        return
    end
 
    datos=filename;
    [A,V,P,t,sensores,G,Mag,tmag]=cadera(datos);
    handles.A=A;
    handles.V=V;
    handles.P=P;
    handles.t=t;
    handles.sensores=sensores;
    handles.G=G;
    handles.Mag=Mag;
    handles.tmag=tmag;
    cd (mainfolder)
    if isnan(str2double(get(handles.nsensor,'string'))) || isnan(str2double(get(handles.tinicial,'string'))) || isnan(str2double(get(handles.tfinal,'string')))
        warndlg(['Si no se modifican, los parámetros por defecto son los siguientes:' sprintf('\n') 'N° Sensor = 1' sprintf('\n') 'Filtro = 0.5' sprintf('\n') 'Tiempo inicial = 0' sprintf('\n') 'Tiempo final = Tiempo máximo del archivo'])
    end
    
    if isnan(str2double(get(handles.nsensor,'string')))
        sensor=1;
    else
        sensor=round(str2double(get(handles.nsensor,'string')));
    end
    if sensor>sensores
        sensor=sensores;
    end
    if sensor<1
        sensor=1;
    end
    handles.sensor=sensor;
    
    if isnan(str2double(get(handles.tinicial,'string')))
        cti=1;
    else
        ti=str2double(get(handles.tinicial,'string'));
        cti=tiempo(ti,t);
    end
    handles.cti=cti;
    
    if isnan(str2double(get(handles.tfinal,'string')))
        ctf=length(t);
    else
        tf=str2double(get(handles.tfinal,'string'));
        ctf=tiempo(tf,t);
    end
    handles.ctf=ctf;
    
    if isnan(str2double(get(handles.filtro,'string')))
        filt=0.5;
    else
        filt=str2double(get(handles.filtro,'string'));
    end
    filt=0.3;
    filtrobajo=0.01;
    [pos,vel,posPlot,quatPlot,gyrX,gyrY,gyrZ,accX,accY,accZ,acc,time,acc_magFilt,stationary,samplePeriod,ang]=script2(filename,sensor,filt,filtrobajo);
    handles.ang=ang;
    handles.filt=filt;
    An=norma(A);
    handles.An=An;
    guidata(gcbo,handles);
end
cd(mainfolder)

% ------------------------------------------------------------
% Read the current directory and sort the names
% ------------------------------------------------------------
function load_listbox(dir_path,handles)
cd (dir_path)
dir_struct = dir(dir_path);
[sorted_names,sorted_index] = sortrows({dir_struct.name}');
handles.file_names = sorted_names;
handles.is_dir = [dir_struct.isdir];
handles.sorted_index = sorted_index;
guidata(handles.figure1,handles)
set(handles.lista_datos,'String',handles.file_names,...
	'Value',1)
set(handles.clickarchivo,'String',['Dirección actual: ' pwd])


% --- Executes during object creation, after setting all properties.
function lista_datos_CreateFcn(hObject, ~, ~)
% hObject    handle to lista_datos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background, change
%       'usewhitebg' to 0 to use default.  See ISPC and COMPUTER.
usewhitebg = 1;
if usewhitebg
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(groot,'defaultUicontrolBackgroundColor'));
end


% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, ~, ~)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Add the current directory to the path, as the pwd might change thru' the
% gui. Remove the directory from the path when gui is closed 
% (See figure1_DeleteFcn)
setappdata(hObject, 'StartPath', pwd);
addpath(pwd);


% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, ~, ~)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Remove the directory added to the path in the figure1_CreateFcn.
if isappdata(hObject, 'StartPath')
    rmpath(getappdata(hObject, 'StartPath'));
end

%%%%%BOTONES%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on button press in animacion. %ANIMACION
function animacion_Callback(~, ~, handles)
% hObject    handle to animacion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%%
%valores por defecto
comp=get(handles.checkbox1,'Value');
filename=handles.filename;
cla(handles.axes1)
sensor=handles.sensor;
sensor2=handles.sensor2;
sensores=handles.sensores;
filt=handles.filt;
if isnan(str2double(get(handles.nsensor,'string')))
    sensor=1;
else
    sensor=round(str2double(get(handles.nsensor,'string')));
end
if sensor>sensores
    sensor=sensores;
end
if sensor<1
    sensor=1;
end
handles.sensor=sensor;
if isnan(str2double(get(handles.sensor2,'string')))
    sensor2=sensores;
else
    sensor2=round(str2double(get(handles.sensor2,'string')));
end
if sensor2>sensores
    sensor2=sensores;
end
if sensor2<1
    sensor2=1;
end
handles.sensor2=sensor2;
if isnan(str2double(get(handles.filtro,'string')))
    filt=0.5;
else
    filt=str2double(get(handles.filtro,'string'));
end
handles.filt=filt;
filtrobajo=0.01;
[pos,vel,posPlot,quatPlot,gyrX,gyrY,gyrZ,accX,accY,accZ,acc,time,acc_magFilt,stationary,samplePeriod,ang]=script2(filename,sensor,filt,filtrobajo);
A=handles.A;
V=handles.V;
P=handles.P;
t=handles.t;
G=handles.G;
Mag=handles.Mag;
tmag=handles.tmag;
cti=handles.cti;
ctf=handles.ctf;
if isnan(str2double(get(handles.tinicial,'string')))
    cti=1;
else
    ti=str2double(get(handles.tinicial,'string'));
    cti=tiempo(ti,time);
end
handles.cti=cti;

if isnan(str2double(get(handles.tfinal,'string')))
    ctf=length(time);
else
    tf=str2double(get(handles.tfinal,'string'));
    ctf=tiempo(tf,time);
end
handles.ctf=ctf;
if comp==0
    t=t(cti:ctf);
    axes(handles.axes1);
    hold on
    plot(t,acc_magFilt(cti:ctf))
    ylabel('Magnitud de la aceleración[g]')
    xlabel('Tiempo[s]')
    xlim([t(cti) t(end)])
    ylim([min(acc_magFilt) max(acc_magFilt)])
    [x y]=ginput(1);
    h1 = text(x,y,'Umbral de paso', ...
        'HorizontalAlignment','center', ...
        'Color', [ 0.3137    0.3137    0.3137], ...
        'FontSize',10);
    filt=y;
    refline(0,y)
    [pos,vel,posPlot,quatPlot,gyrX,gyrY,gyrZ,accX,accY,accZ,acc,time,acc_magFilt,stationary,samplePeriod,ang]=script2(filename,sensor,filt,filtrobajo);
    figure('Position', [9 39 900 600], 'NumberTitle', 'off', 'Name', 'Sensor Data');
    ax(1) = subplot(2,1,1);
    hold on;
    plot(t, gyrX(cti:ctf), 'r');
    plot(t, gyrY(cti:ctf), 'g');
    plot(t, gyrZ(cti:ctf), 'b');
    xlim([t(cti) t(end)]);
    title('Gyroscope');
    xlabel('Time (s)');
    ylabel('Angular velocity (^\circ/s)');
    legend('X', 'Y', 'Z');
    hold off;
    ax(2) = subplot(2,1,2);
    hold on;
    plot(t, accX(cti:ctf), 'r');
    plot(t, accY(cti:ctf), 'g');
    plot(t, accZ(cti:ctf), 'b');
    plot(t, acc_magFilt(cti:ctf), ':k');
    plot(t, stationary(cti:ctf), 'k', 'LineWidth', 2);
    title('Accelerometer');
    xlabel('Time (s)');
    ylabel('Acceleration (g)');
    xlim([t(cti) t(end)]);
    legend('X', 'Y', 'Z', 'Filtered', 'Stationary');
    hold off;
    linkaxes(ax,'x');
    
    %Plot translational accelerations
    figure('Position', [9 39 900 300], 'NumberTitle', 'off', 'Name', 'Accelerations');
    hold on;
    plot(t, acc((cti:ctf),1), 'r');
    plot(t, acc((cti:ctf),2), 'g');
    plot(t, acc((cti:ctf),3), 'b');
    title('Acceleration');
    xlabel('Time (s)');
    ylabel('Acceleration (m/s/s)');
    xlim([t(cti) t(end)]);
    legend('X', 'Y', 'Z');
    hold off;
    
    % Plot translational velocity
    figure('Position', [9 39 900 300], 'NumberTitle', 'off', 'Name', 'Velocity');
    hold on;
    plot(t, vel((cti:ctf),1), 'r');
    plot(t, vel((cti:ctf),2), 'g');
    plot(t, vel((cti:ctf),3), 'b');
    title('Velocity');
    xlabel('Time (s)');
    ylabel('Velocity (m/s)');
    xlim([t(cti) t(end)]);
    legend('X', 'Y', 'Z');
    hold off;
    
    % Plot translational position
    figure('Position', [9 39 900 600], 'NumberTitle', 'off', 'Name', 'Position');
    hold on;
    plot(t, pos((cti:ctf),1), 'r');
    plot(t, pos((cti:ctf),2), 'g');
    plot(t, pos((cti:ctf),3), 'b');
    title('Position');
    xlabel('Time (s)');
    ylabel('Position (m)');
    xlim([t(cti) t(end)]);
    legend('X', 'Y', 'Z');
    hold off;
    cd(handles.mainfolder)
    % Create 6 DOF animation
    SamplePlotFreq = 4;
    Spin = 120;
    posPlot=posPlot(cti:ctf,:);
    quatPlot=quatPlot(cti:ctf,:);  
    SixDofAnimation(posPlot, quatern2rotMat(quatPlot), ...
        'SamplePlotFreq', SamplePlotFreq, 'Trail', 'All', ...
        'Position', [9 39 1280 668], 'View', [(100:(Spin/(length(posPlot)-1)):(100+Spin))', 10*ones(length(posPlot), 1)], ...
        'AxisLength', 0.1, 'ShowArrowHead', false, ...
        'Xlabel', 'X (m)', 'Ylabel', 'Y (m)', 'Zlabel', 'Z (m)','Title',['Sensor ',num2str(sensor)], 'ShowLegend', false, ...
        'CreateAVI', false, 'AVIfileNameEnum', false, 'AVIfps', ((1/samplePeriod) / SamplePlotFreq));
 
else %caso comparación de sensores
    filtrobajo=0.01;
    [pos1,vel1,posPlot1,quatPlot1,gyrX1,gyrY1,gyrZ1,accX1,accY1,accZ1,acc1,time1,acc_magFilt1,stationary1,samplePeriod1,ang1]=script2(filename,sensor,filt,filtrobajo);
    [pos2,vel2,posPlot2,quatPlot2,gyrX2,gyrY2,gyrZ2,accX2,accY2,accZ2,acc2,time2,acc_magFilt2,stationary2,samplePeriod2,ang2]=script2(filename,sensor2,filt,filtrobajo);
    cd(datafolder)
    figure
    hold on
    
    plot(time1(cti:ctf),acc_magFilt1(cti:ctf))
    ylabel('Magnitud de la aceleración[g]')
    xlabel('Tiempo[s]')
    xlim([time1(cti) time1(ctf)])
    ylim([min(acc_magFilt1) max(acc_magFilt1)])
    [x1,y1]=ginput(1);
    h1 = text(x1,y1,'Filtro de paso', ...
        'HorizontalAlignment','center', ...
        'Color', [ 0.3137    0.3137    0.3137], ...
        'FontSize',10);
    filt1=y1;
    refline(0,y1)
    %%

    figure
    hold on
    plot(time2(1:end),acc_magFilt2(1:length(time2)))
    ylabel('Magnitud de la aceleración[g]')
    xlabel('Tiempo[s]')
    xlim([time2(cti) time2(ctf)])
    ylim([min(acc_magFilt2) max(acc_magFilt2)])
    x2=0;
    y2=0;
    [x2 y2]=ginput(1);
    h2 = text(x2,y2,'Filtro de paso', ...
        'HorizontalAlignment','center', ...
        'Color', [ 0.3137    0.3137    0.3137], ...
        'FontSize',10);
    filt2=y2;
    refline(0,y2)
    
    [pos1,vel1,posPlot1,quatPlot1,gyrX1,gyrY1,gyrZ1,accX1,accY1,accZ1,acc1,time1,acc_magFilt1,stationary1,samplePeriod1]=script2(filename,sensor,filt1,filtrobajo);
    [pos2,vel2,posPlot2,quatPlot2,gyrX2,gyrY2,gyrZ2,accX2,accY2,accZ2,acc2,time2,acc_magFilt2,stationary2,samplePeriod2]=script2(filename,sensor2,filt2,filtrobajo);
    %%
    figure('Position', [9 39 900 600], 'NumberTitle', 'off', 'Name', 'Sensor Data');
ax(1) = subplot(2,1,1);
    hold on;
    plot(time1, accX1, 'r');
    plot(time1, accY1, 'g');
    plot(time1, accZ1, 'b');
    plot(time1, acc_magFilt1, ':k');
    plot(time1, stationary1, 'k', 'LineWidth', 2);
    title('Sensor 1');
    xlabel('Time (s)');
    ylabel('Aceleración (g)');
    legend('X', 'Y', 'Z', 'Filtered', 'Stationary');
    hold off;
ax(2) = subplot(2,1,2);
    hold on;
    plot(time2, accX2, 'r');
    plot(time2, accY2, 'g');
    plot(time2, accZ2, 'b');
    plot(time2, acc_magFilt2, ':k');
    plot(time2, stationary2, 'k', 'LineWidth', 2);
    title('Sensor 2');
    xlabel('Time (s)');
    ylabel('Aceleración (g)');
    legend('X', 'Y', 'Z', 'Filtered', 'Stationary');
    hold off;
linkaxes(ax,'x');

figure('Position', [9 39 900 600], 'NumberTitle', 'off', 'Name', 'Sensor Data');
ax(1) = subplot(2,1,1);
    hold on;
    plot(time1, gyrX1, 'r');
    plot(time1, gyrY1, 'g');
    plot(time1, gyrZ1, 'b');
    title('Gyroscope');
    xlabel('Time (s)');
    ylabel('Angular velocity (^\circ/s)');
    legend('X', 'Y', 'Z');
    hold off;
ax(2) = subplot(2,1,2);
    hold on;
    plot(time2, gyrX2, 'r');
    plot(time2, gyrY2, 'g');
    plot(time2, gyrZ2, 'b');
    title('Gyroscope');
    xlabel('Time (s)');
    ylabel('Angular velocity (^\circ/s)');
    legend('X', 'Y', 'Z');
    hold off;
linkaxes(ax,'x');
    %%
    %%
    SamplePlotFreq = 4;
    Spin = 120;
    [pos1,vel1,posPlot1,quatPlot1,gyrX1,gyrY1,gyrZ1,accX1,accY1,accZ1,acc1,time1,acc_magFilt1,stationary1,samplePeriod1]=script2(filename,sensor,filt1,filtrobajo);
    [pos2,vel2,posPlot2,quatPlot2,gyrX2,gyrY2,gyrZ2,accX2,accY2,accZ2,acc2,time2,acc_magFilt2,stationary2,samplePeriod2]=script2(filename,sensor2,filt2,filtrobajo);
    %%
    comp
        posPlot1=posPlot1(cti:ctf,:);
        quatPlot1=quatPlot1(cti:ctf,:);
        SixDofAnimation(posPlot1, quatern2rotMat(quatPlot1), ...
        'SamplePlotFreq', SamplePlotFreq, 'Trail', 'All', ...
        'Position', [9 39 1280 768], 'View', [(100:(Spin/(length(posPlot1)-1)):(100+Spin))', 10*ones(length(posPlot1), 1)], ...
        'AxisLength', 0.1, 'ShowArrowHead', false, ...
        'Xlabel', 'X (m)', 'Ylabel', 'Y (m)', 'Zlabel', 'Z (m)','Title',['Sensor ',num2str(sensor)], 'ShowLegend', false, ...
        'CreateAVI', false, 'AVIfileNameEnum', false, 'AVIfps', ((1/samplePeriod1) / SamplePlotFreq));
    %%
        posPlot2=posPlot2(cti:ctf,:);
        quatPlot2=quatPlot2(cti:ctf,:);
        SixDofAnimation(posPlot2, quatern2rotMat(quatPlot2), ...
        'SamplePlotFreq', SamplePlotFreq, 'Trail', 'All', ...
        'Position', [9 39 1280 768], 'View', [(100:(Spin/(length(posPlot2)-1)):(100+Spin))', 10*ones(length(posPlot2), 1)], ...
        'AxisLength', 0.1, 'ShowArrowHead', false, ...
        'Xlabel', 'X (m)', 'Ylabel', 'Y (m)', 'Zlabel', 'Z (m)','Title',['Sensor ',num2str(sensor2)], 'ShowLegend', false, ...
        'CreateAVI', false, 'AVIfileNameEnum', false, 'AVIfps', ((1/samplePeriod2) / SamplePlotFreq));
end          
% --- Executes on button press in grafico_dimensional.
function grafico_dimensional_Callback(~, ~, handles) %grafico dimensional
% hObject    handle to grafico_dimensional (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%%
%valores por defecto
A=handles.A;
V=handles.V;
P=handles.P;
t=handles.t;
sensores=handles.sensores;
G=handles.G;
Mag=handles.Mag;
tmag=handles.tmag;
sensor=handles.sensor;
cti=handles.cti;
ctf=handles.ctf;
filt=handles.filt;

if isnan(str2double(get(handles.nsensor,'string')))
    sensor=1;
else
    sensor=round(str2double(get(handles.nsensor,'string')));
end
if sensor>sensores
    sensor=sensores;
end
if sensor<1
    sensor=1;
end
handles.sensor=sensor;

if isnan(str2double(get(handles.tinicial,'string')))
    cti=1;
else
    ti=str2double(get(handles.tinicial,'string'));
    cti=tiempo(ti,t);
end
handles.cti=cti;
[l p]=cortarmatrices(A,sensor);
if isnan(str2double(get(handles.tfinal,'string')))
    ctf=p;
else
    tf=str2double(get(handles.tfinal,'string'));
    ctf=tiempo(tf,t);
end
handles.ctf=ctf;

if isnan(str2double(get(handles.filtro,'string')))
    filt=0.5;
else
    filt=str2double(get(handles.filtro,'string'));
end
handles.filt=filt;
%%

cla(handles.axes1)

axes(handles.axes1);
plot(t(cti:ctf),A(1,cti:ctf,sensor),'r')
hold on
plot(t(cti:ctf),A(2,cti:ctf,sensor),'g')
plot(t(cti:ctf),A(3,cti:ctf,sensor),'b')
legend('X','Y','Z')
xlim([t(cti) t(ctf)])
ylim([min(min(min(A(:,cti:ctf,sensor))))-0.5 max(max(max(A(:,cti:ctf,sensor))))+0.5])
title(['Sensor ',num2str(sensor)])

xlabel('Tiempo[s]')
ylabel('Aceleración[g]')
[x,~]=ginput(2);
cti=tiempo(min(x),t);
ctf=tiempo(max(x),t);
hold off
plot(t(cti:ctf),A(1,cti:ctf,sensor),'r')
hold on
plot(t(cti:ctf),A(2,cti:ctf,sensor),'g')
plot(t(cti:ctf),A(3,cti:ctf,sensor),'b')
legend('X','Y','Z')
xlim([t(cti) t(ctf)])
ylim([min(min(min(A(:,cti:ctf,sensor))))-0.5 max(max(max(A(:,cti:ctf,sensor))))+0.5])
title(['Sensor ',num2str(sensor)])
xlabel('Tiempo[s]')
ylabel('Aceleración[g]')

%%
%plotear velocidad y posición
% set(gca,'XTick',[])
% axes(handles.axes2);
% plot(t(cti:ctf),V(1,cti:ctf,sensor).*9.81,'r')
% hold on
% plot(t(cti:ctf),V(2,cti:ctf,sensor).*9.81,'b')
% plot(t(cti:ctf),V(3,cti:ctf,sensor).*9.81,'g')
% %xlabel('Tiempo[s]')
% ylabel('Velocidad[m/s]')
% set(gca,'XTick',[])
% %legend('Velocidad x','Velocidad y','Velocidad z')
% axes(handles.axes3);
% plot(t(cti:ctf),P(1,cti:ctf,sensor).*9.81,'r')
% hold on
% plot(t(cti:ctf),P(2,cti:ctf,sensor).*9.81,'b')
% plot(t(cti:ctf),P(3,cti:ctf,sensor).*9.81,'g')
% xlabel('Tiempo[s]')
% ylabel('Posición[m]')
% %legend('Posición x','Posición y','Posición z')
%%
guidata(gcbo,handles);


% --- Executes on button press in grafico_magnitud.
function grafico_magnitud_Callback(~, ~, handles) %magnitud
% hObject    handle to grafico_magnitud (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%%
%valores por defecto
A=handles.A;
V=handles.V;
P=handles.P;
t=handles.t;
sensores=handles.sensores;
G=handles.G;
Mag=handles.Mag;
tmag=handles.tmag;
sensor=handles.sensor;
sensor2=handles.sensor2;
cti=handles.cti;
ctf=handles.ctf;
filt=handles.filt;
if isnan(str2double(get(handles.nsensor,'string')))
    sensor=1;
else
    sensor=round(str2double(get(handles.nsensor,'string')));
end
if sensor>sensores
    sensor=sensores;
end
if sensor<1
    sensor=1;
end
handles.sensor=sensor;

if isnan(str2double(get(handles.tinicial,'string')))
    cti=1;
else
    ti=str2double(get(handles.tinicial,'string'));
    cti=tiempo(ti,t);
end
handles.cti=cti;

[l p]=cortarmatrices(A,sensor);

if isnan(str2double(get(handles.tfinal,'string')))
    ctf=p;
else
    tf=str2double(get(handles.tfinal,'string'));
    ctf=tiempo(tf,t);
end
handles.ctf=ctf;

if isnan(str2double(get(handles.filtro,'string')))
    filt=0.5;
else
    filt=str2double(get(handles.filtro,'string'));
end
handles.filt=filt;
%%

cla(handles.axes1)

An=norma(A);
Vn=norma(V);
Pn=norma(P);

axes(handles.axes1);
legend(handles.axes1,'hide')
plot(t(cti:ctf),An(cti:ctf,sensor),'r')
title(['Sensor ',num2str(sensor)])
xlim([t(cti) t(ctf)])
ylim([min(min(min(An(cti:ctf,sensor))))-0.5 max(max(max(An(cti:ctf,sensor))))+0.5])
xlabel('Tiempo [s]')
ylabel('Aceleración[g]')
[x , ~]=ginput(2);
cti=tiempo(min(x),t);
ctf=tiempo(max(x),t);
hold off
legend(handles.axes1,'hide')
plot(t(cti:ctf),An(cti:ctf,sensor),'r')
xlim([t(cti) t(ctf)])
ylim([min(min(min(An(cti:ctf,sensor))))-0.5 max(max(max(An(cti:ctf,sensor))))+0.5])
title(['Sensor ',num2str(sensor)])
xlabel('Tiempo [s]')
ylabel('Aceleración[g]')
handles.An=An;

%%
%plotear velocidad y posición
% set(gca,'XTick',[])
% axes(handles.axes2);
% plot(t1(cti1:ctf1),Vn(cti1:ctf1,sensor).*9.81,'g')
% ylabel('Velocidad[m/s]')
% set(gca,'XTick',[])
% axes(handles.axes3);
% plot(t1(cti1:ctf1),Pn(cti1:ctf1,sensor).*9.81,'b')
% xlabel('Tiempo[s]')
% ylabel('Posición[m]')
%%
guidata(gcbo,handles);



% --- Executes on button press in giroscopio.
function giroscopio_Callback(~, ~, handles) %giroscopio
% hObject    handle to giroscopio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%%
%valores por defecto
A=handles.A;
V=handles.V;
P=handles.P;
t=handles.t;
sensores=handles.sensores;
G=handles.G;
Mag=handles.Mag;
tmag=handles.tmag;
sensor=handles.sensor;
cti=handles.cti;
ctf=handles.ctf;
filt=handles.filt;
if isnan(str2double(get(handles.nsensor,'string')))
    sensor=1;
else
    sensor=round(str2double(get(handles.nsensor,'string')));
end
if sensor>sensores
    sensor=sensores;
end
if sensor<1
    sensor=1;
end
handles.sensor=sensor;

if isnan(str2double(get(handles.tinicial,'string')))
    cti=1;
else
    ti=str2double(get(handles.tinicial,'string'));
    cti=tiempo(ti,t);
end
handles.cti=cti;

[l p]=cortarmatrices(A,sensor);

if isnan(str2double(get(handles.tfinal,'string')))
    ctf=p;
else
    tf=str2double(get(handles.tfinal,'string'));
    ctf=tiempo(tf,t);
end
handles.ctf=ctf;

if isnan(str2double(get(handles.filtro,'string')))
    filt=0.5;
else
    filt=str2double(get(handles.filtro,'string'));
end
handles.filt=filt;
%%

cla(handles.axes1)

axes(handles.axes1);
plot(t(cti:ctf),G(1,cti:ctf,sensor),'r')
hold on
plot(t(cti:ctf),G(2,cti:ctf,sensor),'g')
plot(t(cti:ctf),G(3,cti:ctf,sensor),'b')
legend('X','Y','Z')
title(['Sensor ',num2str(sensor)])
xlabel('Tiempo[s]')
ylabel('Vel. Angular[°/s]')
xlim([t(cti) t(ctf)])
ylim([min(min(min(G(:,cti:ctf,sensor))))-20 max(max(max(G(:,cti:ctf,sensor))))+20])
[x , ~]=ginput(2);
hold off
cti=tiempo(min(x),t);
ctf=tiempo(max(x),t);
plot(t(cti:ctf),G(1,cti:ctf,sensor),'r')
hold on
plot(t(cti:ctf),G(2,cti:ctf,sensor),'b')
plot(t(cti:ctf),G(3,cti:ctf,sensor),'g')
legend('X','Y','Z')
title(['Sensor ',num2str(sensor)])
xlabel('Tiempo[s]')
ylabel('Vel. Angular[°/s]')
xlim([t(cti) t(ctf)])
ylim([min(min(min(G(:,cti:ctf,sensor))))-20 max(max(max(G(:,cti:ctf,sensor))))+20])

guidata(gcbo,handles);


% --- Executes on button press in magnetometro.
function magnetometro_Callback(~, ~, handles) %magnetómetro
% hObject    handle to magnetometro (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

A=handles.A;
V=handles.V;
P=handles.P;
t=handles.t;
sensores=handles.sensores;
G=handles.G;
Mag=handles.Mag;
tmag=handles.tmag;
sensor=handles.sensor;
sensor2=handles.sensor2;
cti=handles.cti;
ctf=handles.ctf;
filt=handles.filt;
if isnan(str2double(get(handles.nsensor,'string')))
    sensor=1;
else
    sensor=round(str2double(get(handles.nsensor,'string')));
end
if isnan(str2double(get(handles.sensor2,'string')))
    sensor2=1;
else
    sensor2=round(str2double(get(handles.sensor2,'string')));
end
if sensor2>sensores
    sensor2=sensores;
end
if sensor2<1
    sensor2=1;
end
if sensor>sensores
    sensor=sensores;
end
if sensor<1
    sensor=1;
end
handles.sensor=sensor;
handles.sensor2=sensor2;
if isnan(str2double(get(handles.tinicial,'string')))
    cti=1;
else
    ti=str2double(get(handles.tinicial,'string'));
    cti=tiempo(ti,tmag);
end
handles.cti=cti;

[l p]=cortarmatrices(Mag,sensor);

if isnan(str2double(get(handles.tfinal,'string')))
    ctf=p;
else
    tf=str2double(get(handles.tfinal,'string'));
    ctf=tiempo(tf,tmag);
end
handles.ctf=ctf;

if isnan(str2double(get(handles.filtro,'string')))
    filt=0.5;
else
    filt=str2double(get(handles.filtro,'string'));
end
handles.filt=filt;
%%

cla(handles.axes1)


axes(handles.axes1);
plot(tmag(cti:ctf),Mag(1,cti:ctf,sensor),'r')
hold on
plot(tmag(cti:ctf),Mag(2,cti:ctf,sensor),'g')
plot(tmag(cti:ctf),Mag(3,cti:ctf,sensor),'b')
legend('X','Y','Z')
title(['Sensor ',num2str(sensor)])
xlabel('Tiempo[s]')
ylabel('miliTesla')
xlim([tmag(cti) tmag(ctf)])
ylim([min(min(min(Mag(:,cti:ctf,sensor)))) max(max(max(Mag(:,cti:ctf,sensor))))])
[x , ~]=ginput(2);
cti=tiempo(min(x),tmag);
ctf=tiempo(max(x),tmag);
hold off
plot(tmag(cti:ctf),Mag(1,cti:ctf,sensor),'r')
hold on
plot(tmag(cti:ctf),Mag(2,cti:ctf,sensor),'b')
plot(tmag(cti:ctf),Mag(3,cti:ctf,sensor),'g')
legend('X','Y','Z')
title(['Sensor ',num2str(sensor)])
xlabel('Tiempo[s]')
ylabel('miliTesla')
xlim([tmag(cti) tmag(ctf)])
ylim([min(min(min(Mag(:,cti:ctf,sensor)))) max(max(max(Mag(:,cti:ctf,sensor))))])

guidata(gcbo,handles);



% --- Executes on button press in angulos.
function angulos_Callback(hObject, eventdata, handles)
% hObject    handle to angulos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
comp=get(handles.checkbox1,'Value');
filename=handles.filename;
cla(handles.axes1)
sensor=handles.sensor;
sensor2=handles.sensor2;
sensores=handles.sensores;
filt=handles.filt;
A=handles.A;
V=handles.V;
P=handles.P;
t=handles.t;
sensores=handles.sensores;
G=handles.G;
Mag=handles.Mag;
tmag=handles.tmag;
sensor=handles.sensor;
sensor2=handles.sensor2;
cti=handles.cti;
ctf=handles.ctf;
filt=handles.filt;
if isnan(str2double(get(handles.nsensor,'string')))
    sensor=1;
else
    sensor=round(str2double(get(handles.nsensor,'string')));
end
if sensor>sensores
    sensor=sensores;
end
if sensor<1
    sensor=1;
end
handles.sensor=sensor;

if isnan(str2double(get(handles.tinicial,'string')))
    cti=1;
else
    ti=str2double(get(handles.tinicial,'string'));
    cti=tiempo(ti,t);
end
handles.cti=cti;

[l p]=cortarmatrices(A,sensor);

if isnan(str2double(get(handles.tfinal,'string')))
    ctf=p;
else
    tf=str2double(get(handles.tfinal,'string'));
    ctf=tiempo(tf,t);
end
handles.ctf=ctf;

if isnan(str2double(get(handles.filtro,'string')))
    filt=0.5;
else
    filt=str2double(get(handles.filtro,'string'));
end
handles.filt=filt;
%%

cla(handles.axes1)

An=norma(A);
Vn=norma(V);
Pn=norma(P);

axes(handles.axes1);
legend(handles.axes1,'hide')
plot(t(cti:ctf),An(cti:ctf,sensor),'r')
title(['Sensor ',num2str(sensor)])
xlim([t(cti) t(ctf)])
ylim([min(min(min(An(cti:ctf,sensor))))-0.5 max(max(max(An(cti:ctf,sensor))))+0.5])
xlabel('Tiempo [s]')
ylabel('Aceleración[g]')
hold off
legend(handles.axes1,'hide')
plot(t(cti:ctf),An(cti:ctf,sensor),'r')
xlim([t(cti) t(ctf)])
ylim([min(min(min(An(cti:ctf,sensor))))-0.5 max(max(max(An(cti:ctf,sensor))))+0.5])
title(['Sensor ',num2str(sensor)])
xlabel('Tiempo [s]')
ylabel('Aceleración[g]')
[Tiempoangulo,~]=getpts(handles.axes1);
filtrobajo=0.01;
filt=0.3;
[pos,vel,posPlot,quatPlot,gyrX,gyrY,gyrZ,accX,accY,accZ,acc,time,acc_magFilt,stationary,samplePeriod,ang]=script2(filename,sensor,filt,filtrobajo);
tiempoangulofinal=tiempo(Tiempoangulo(end),t);
x=num2str(ang(tiempoangulofinal,1));
y=num2str(ang(tiempoangulofinal,2));
z=num2str(ang(tiempoangulofinal,3));
msgbox(['giro en X ' x sprintf('\n') 'giro en Y ' y sprintf('\n') 'giro en Z ' z])
handles.An=An;
handles.ang=ang;
guidata(gcbo,handles);


% --- Executes on button press in reset.
function reset_Callback(~, ~, handles) %reset
% hObject    handle to reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

cla(handles.axes1)
comp=get(handles.checkbox1,'Value');

if comp==0
    close figure 1
    close figure 2
    close figure 3
    close figure 4
    close figure 5
else
    close figure 1
    close figure 2
    close figure 3
    close figure 4
    close figure 5
    close figure 6
end

% --- Executes on button press in info.
function info_Callback(hObject, eventdata, handles)
% hObject    handle to info (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cd(handles.mainfolder)
h = msgbox(['1) Se deben agregar las carpetas con archivos al ''path'' de Matlab.' ...
    sprintf('\n') 'HOME->Set Path->Add with subfolders->Seleccionar carpetas deseadas' ...
    sprintf('\n') sprintf('\n') '2) Parámetros por defecto:' sprintf('\n') ...
    'N° Sensor = 1' sprintf('\n') 'Filtro = 0.5' sprintf('\n') 'Tiempo inicial = 0' ...
    sprintf('\n')  'Tiempo final = Tiempo máximo del archivo' sprintf('\n') sprintf('\n')...
    '3) Se recomienda tener una carpeta con todos los archivos .m y dentro de' ...
    'esta una carpeta ''Datasets'' con los archivos .csv'],'Información','custom.jpg');

%%%EDIT TEXTS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function tinicial_Callback(hObject, ~, handles)
% hObject    handle to tinicial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tinicial as text
%        str2double(get(hObject,'String')) returns contents of tinicial as a double
handles.tinicial=str2double(get(hObject,'String'));


% --- Executes during object creation, after setting all properties.
function tinicial_CreateFcn(hObject, ~, ~)
% hObject    handle to tinicial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tfinal_Callback(hObject, ~, handles)
% hObject    handle to tfinal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tfinal as text
%        str2double(get(hObject,'String')) returns contents of tfinal as a double
handles.tfinal=str2double(get(hObject,'String'));


% --- Executes during object creation, after setting all properties.
function tfinal_CreateFcn(hObject, ~, ~)
% hObject    handle to tfinal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function nsensor_Callback(hObject, ~, handles)
% % hObject    handle to nsensor (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% % Hints: get(hObject,'String') returns contents of nsensor as text
% %        str2double(get(hObject,'String')) returns contents of nsensor as a double
handles.nsensor=str2double(get(hObject,'String'));


% --- Executes during object creation, after setting all properties.
function nsensor_CreateFcn(hObject, ~, ~)
% hObject    handle to nsensor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function filtro_Callback(~, ~, ~)
% hObject    handle to filtro (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of filtro as text
%        str2double(get(hObject,'String')) returns contents of filtro as a double


% --- Executes during object creation, after setting all properties.
function filtro_CreateFcn(hObject, ~, ~)
% hObject    handle to filtro (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1



function sensor2_Callback(hObject, eventdata, handles)
% hObject    handle to sensor2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sensor2 as text
%        str2double(get(hObject,'String')) returns contents of sensor2 as a double
handles.sensor2=str2double(get(hObject,'String'));


% --- Executes during object creation, after setting all properties.
function sensor2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sensor2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function directorio_principal_Callback(hObject, eventdata, handles)
% hObject    handle to directorio_principal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of directorio_principal as text
%        str2double(get(hObject,'String')) returns contents of directorio_principal as a double
handles.mainfolder=get(hObject,'String');
cd(handles.mainfolder)
% --- Executes during object creation, after setting all properties.
function directorio_principal_CreateFcn(hObject, eventdata, handles)
% hObject    handle to directorio_principal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function directorio_datos_Callback(hObject, eventdata, handles)
% hObject    handle to directorio_datos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of directorio_datos as text
%        str2double(get(hObject,'String')) returns contents of directorio_datos as a double
cd(get(hObject,'String'))
handles.datafolder=get(hObject,'String');
if nargin == 3
    initial_dir = pwd;
elseif nargin > 4
    if strcmpi(varargin{1},'dir')
        if exist(varargin{2},'dir')
            initial_dir = varargin{2};
        else
            errordlg('Input argument must be a valid directory','Input Argument Error!')
            return
        end
    else
        errordlg('Unrecognized input argument','Input Argument Error!');
        return;
    end
end
load_listbox(initial_dir,handles)


% --- Executes during object creation, after setting all properties.
function directorio_datos_CreateFcn(hObject, eventdata, handles)
% hObject    handle to directorio_datos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end






% --- Executes on button press in exportar.
function exportar_Callback(hObject, eventdata, handles)
% hObject    handle to exportar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
An=handles.An;
ang=handles.ang;
t=handles.t;
H=get(handles.dataexport,'String');
sensor=handles.sensor;
[a b]=size(t);
L=zeros(a-length(ang(:,1)),3);
Ang=[ang;L];
titulo={'Tiempo','Aceleracion','AnguloX','AnguloY','AnguloZ'};
S=[t An(:,sensor) Ang];
sTable=array2table(S,'VariableNames',titulo);
writetable(sTable,[H '.xls']);



function dataexport_Callback(hObject, eventdata, handles)
% hObject    handle to dataexport (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dataexport as text
%        str2double(get(hObject,'String')) returns contents of dataexport as a double
handles.export=get(hObject,'String');

% --- Executes during object creation, after setting all properties.
function dataexport_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dataexport (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
