function varargout = Hydraulik_YouTube_02_GUI(varargin)
% HYDRAULIK_YOUTUBE_02_GUI MATLAB code for Hydraulik_YouTube_02_GUI.fig
%      HYDRAULIK_YOUTUBE_02_GUI, by itself, creates a new HYDRAULIK_YOUTUBE_02_GUI or raises the existing
%      singleton*.
%
%      H = HYDRAULIK_YOUTUBE_02_GUI returns the handle to a new HYDRAULIK_YOUTUBE_02_GUI or the handle to
%      the existing singleton*.
%
%      HYDRAULIK_YOUTUBE_02_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in HYDRAULIK_YOUTUBE_02_GUI.M with the given input arguments.
%
%      HYDRAULIK_YOUTUBE_02_GUI('Property','Value',...) creates a new HYDRAULIK_YOUTUBE_02_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Hydraulik_YouTube_02_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Hydraulik_YouTube_02_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Hydraulik_YouTube_02_GUI

% Last Modified by GUIDE v2.5 03-Aug-2020 02:18:27

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Hydraulik_YouTube_02_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @Hydraulik_YouTube_02_GUI_OutputFcn, ...
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


% --- Executes just before Hydraulik_YouTube_02_GUI is made visible.
function Hydraulik_YouTube_02_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Hydraulik_YouTube_02_GUI (see VARARGIN)

global gstop
gstop = 0;

global w1
global w2
global w3
w1 = round(get(handles.slider_w1,'Value'),1);
w2 = round(get(handles.slider_w2,'Value'),1);
w3 = round(get(handles.slider_w3,'Value'),1);
set(handles.text_w1,'String',['w1 = ',num2str(w1)])
set(handles.text_w2,'String',['w2 = ',num2str(w2)])
set(handles.text_w3,'String',['w3 = ',num2str(w3)])

% Choose default command line output for Hydraulik_YouTube_02_GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Hydraulik_YouTube_02_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Hydraulik_YouTube_02_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function slider_w1_Callback(hObject, eventdata, handles)
% hObject    handle to slider_w1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
global w1
w1 = round(get(handles.slider_w1,'Value'),1);
set(handles.text_w1,'String',['w1 = ',num2str(w1)])


% --- Executes during object creation, after setting all properties.
function slider_w1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_w1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider_w2_Callback(hObject, eventdata, handles)
% hObject    handle to slider_w2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
global w2
w2 = round(get(handles.slider_w2,'Value'),1);
set(handles.text_w2,'String',['w2 = ',num2str(w2)])


% --- Executes during object creation, after setting all properties.
function slider_w2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_w2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider_w3_Callback(hObject, eventdata, handles)
% hObject    handle to slider_w3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
global w3
w3 = round(get(handles.slider_w3,'Value'),1);
set(handles.text_w3,'String',['w3 = ',num2str(w3)])


% --- Executes during object creation, after setting all properties.
function slider_w3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_w3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function edit_time_Callback(hObject, eventdata, handles)
% hObject    handle to edit_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_time as text
%        str2double(get(hObject,'String')) returns contents of edit_time as a double


% --- Executes during object creation, after setting all properties.
function edit_time_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_d1_Callback(hObject, eventdata, handles)
% hObject    handle to edit_d1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_d1 as text
%        str2double(get(hObject,'String')) returns contents of edit_d1 as a double


% --- Executes during object creation, after setting all properties.
function edit_d1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_d1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_d2_Callback(hObject, eventdata, handles)
% hObject    handle to edit_d2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_d2 as text
%        str2double(get(hObject,'String')) returns contents of edit_d2 as a double


% --- Executes during object creation, after setting all properties.
function edit_d2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_d2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_d3_Callback(hObject, eventdata, handles)
% hObject    handle to edit_d3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_d3 as text
%        str2double(get(hObject,'String')) returns contents of edit_d3 as a double


% --- Executes during object creation, after setting all properties.
function edit_d3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_d3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pb_start.
function pb_start_Callback(hObject, eventdata, handles)
% hObject    handle to pb_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global w1
global w2
global w3
global gstop
d0 = 5;
d1 = str2double(get(handles.edit_d1,'String'));
d2 = str2double(get(handles.edit_d2,'String'));
d3 = str2double(get(handles.edit_d3,'String'));
A0 = d0^2*pi()/4;
A1 = d1^2*pi()/4;
A2 = d2^2*pi()/4;
A3 = d3^2*pi()/4;
h0 = 4;
h3 = 3;

h = zeros(1)+2;
t = zeros(1);
te = str2double(get(handles.edit_time,'String'));
dt = 0.1;
n = 2;

n_ueber = 1;
h_ueber = [];
t_ueber = [];

ax1 = handles.axes1;
ax2 = handles.axes2;

while t <= te & gstop == 0
    if h(n-1) < h3
        h(n) = h(n-1)+(A1*w1-A2*w2)/A0*dt;
    elseif (h(n-1) >= h3 && h(n-1) < h0) || (h(n-1) >= h0 && A1*w1 < A2*w2+A3*w3)
        h(n) = h(n-1)+(A1*w1-A2*w2-A3*w3)/A0*dt;
    elseif h(n-1) >= h0
        h(n) = h0;
    else
        h(n) = 0;
    end
    t(n) = t(n-1)+dt;
    
    L = d0*0.1;
    Lt = L*0.5;
    plot(ax1,[-d0/2-L d0/2+L],[0 0],'k','LineWidth',3)
    hold(ax1,'on')
    plot(ax1,[-d0/2 -d0/2],[d1 h0],'k','LineWidth',3)
    plot(ax1,[-d0/2-L -d0/2],[d1 d1],'k','LineWidth',3)
    plot(ax1,[d0/2 d0/2],[d2 h3],'k','LineWidth',3)
    plot(ax1,[d0/2 d0/2+L],[d2 d2],'k','LineWidth',3)
    plot(ax1,[d0/2 d0/2],[h3+d3 h0],'k','LineWidth',3)
    plot(ax1,[d0/2 d0/2+L],[h3 h3],'k','LineWidth',3)
    plot(ax1,[d0/2 d0/2+L],[h3+d3 h3+d3],'k','LineWidth',3)
    if h(end) >= h0
        Col = '#D95319';
    else
        Col = '#0072BD';
    end
    plot(ax1,[-d0/2 d0/2],[h(end) h(end)],'Color',Col,'LineWidth',5)
    text(ax1,-d0/2-Lt,d1/2,'1 \rightarrow')
    text(ax1,d0/2-Lt,d2/2,'\rightarrow 2')
    text(ax1,d0/2-Lt,h3+d3/2,'\rightarrow 3')
    ax1.YLim = ([0 h0*1.2]);
    title(ax1,{'Animation'},'FontSize',16,'FontWeight','normal')
    hold(ax1,'off')

    p = plot(ax2,t,h,'k');
    hold(ax2,'on')
    if (length(h) > 1 && h(end) >= h3 && h(end-1) < h3) || (length(h) > 1 && h(end) <= h3 && h(end-1) > h3)
        h_ueber(n_ueber) = h(end);
        t_ueber(n_ueber) = t(end);
        n_ueber = n_ueber+1;
    end
    for m = 1:1:length(h_ueber)
        plot(ax2,[t_ueber(m) t_ueber(m)],[0 h_ueber(m)],'--k')
        plot(ax2,t_ueber(m),h_ueber(m),'k','Marker','.','MarkerSize',15)    
    end
    ax2.XLim = ([0 te]);
    ax2.YLim = ([0 h0*1.2]);
    title(ax2,{'h - t'},'FontSize',16,'FontWeight','normal')
    legend(p,'h(t)','location','NorthEast')
    xlabel('t [s]');
    ylabel('h(t) [m]');
    grid on
    grid minor
    hold(ax2,'off')

    drawnow
   
    n = n+1;
end


% --- Executes on button press in pb_stop.
function pb_stop_Callback(hObject, eventdata, handles)
% hObject    handle to pb_stop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global gstop
gstop = 1;

% --- Executes on button press in pb_close.
function pb_close_Callback(hObject, eventdata, handles)
% hObject    handle to pb_close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pb_stop_Callback()
close

% --- Executes during object creation, after setting all properties.
function text_w1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_w1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function text_w2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_w2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function text_w3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_w3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
