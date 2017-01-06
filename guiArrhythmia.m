function varargout = guiArrhythmia(varargin)
% GUIARRHYTHMIA MATLAB code for guiArrhythmia.fig
%      GUIARRHYTHMIA, by itself, creates a new GUIARRHYTHMIA or raises the existing
%      singleton*.
%
%      H = GUIARRHYTHMIA returns the handle to a new GUIARRHYTHMIA or the handle to
%      the existing singleton*.
%
%      GUIARRHYTHMIA('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUIARRHYTHMIA.M with the given input arguments.
%
%      GUIARRHYTHMIA('Property','Value',...) creates a new GUIARRHYTHMIA or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before guiArrhythmia_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to guiArrhythmia_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help guiArrhythmia

% Last Modified by GUIDE v2.5 06-Jan-2017 21:37:20

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @guiArrhythmia_OpeningFcn, ...
                   'gui_OutputFcn',  @guiArrhythmia_OutputFcn, ...
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


% --- Executes just before guiArrhythmia is made visible.
function guiArrhythmia_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to guiArrhythmia (see VARARGIN)

% Choose default command line output for guiArrhythmia
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes guiArrhythmia wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = guiArrhythmia_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in chooseFile.
function chooseFile_Callback(hObject, eventdata, handles)
% hObject    handle to chooseFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[FileName,PathName,FilterIndex] = uigetfile('*.txt','Choose a file...');
data = dlmread([PathName, FileName],',');
if size(data, 2) == 1
    data = dlmread([PathName, FileName],' ');
end

handles.time = data(1:end,1);
handles.ecg = data(:,2);

handles.fs = 1/mean(diff(handles.time));

% Due to poor resolution of the time samples it is best to calculate a new
% time vector using accurate sampling frequency.
handles.timeNew = (0:length(handles.ecg)-1)'./handles.fs;

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in analyze.
function analyze_Callback(hObject, eventdata, handles)
% hObject    handle to analyze (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[ handles.ecgPreFiltered, handles.highpassed9HzEcg, handles.highpassed3HzEcg, handles.lowpassed5HzEcg ] = ...
    prefilter( handles.ecg, handles.fs );

figure(1)
plot(handles.timeNew,handles.ecgPreFiltered,'k')
hold on
plot(handles.timeNew,handles.highpassed9HzEcg,'b')
plot(handles.timeNew,handles.highpassed3HzEcg,'m')
plot(handles.timeNew,handles.lowpassed5HzEcg,'g')
hold off
grid on


[ handles.rIndeces, handles.rAmplitudes, handles.adaptiveThreshold ] = ...
    detectRPeaks( handles.timeNew, handles.ecgPreFiltered, handles.highpassed9HzEcg, handles.highpassed3HzEcg, handles.lowpassed5HzEcg, handles.fs );

[ handles.rrIntervalsInMs, handles.rPeakTimeStamps ] = calculateRrIntervals( handles.timeNew, handles.rIndeces );

Number_of_approved_R_peaks = length(handles.rrIntervalsInMs)+1;


[ handles.rrCategory ] = classifyRrIntervals( handles.rrIntervalsInMs );

number_of_normal_beats = size(handles.rrCategory(handles.rrCategory == 1), 1);
number_of_premature_ventricular_contractions = size(handles.rrCategory(handles.rrCategory == 2), 1);
number_of_ventricular_flutter_beats = size(handles.rrCategory(handles.rrCategory == 3), 1);
number_of_second_degree_heart_block_beats = size(handles.rrCategory(handles.rrCategory == 4), 1);

sprintf('Number of approved R-peaks: %d\nNumber of normal beats: %d\nNumber of premature ventricular contractions: %d\nNumber of ventricular flutter beats: %d\nNumber of 2nd degree heart block beats: %d',...
    Number_of_approved_R_peaks,number_of_normal_beats,number_of_premature_ventricular_contractions,number_of_ventricular_flutter_beats,number_of_second_degree_heart_block_beats)

set(handles.apprRPeaks,'string',num2str(Number_of_approved_R_peaks))
set(handles.nOfNormal,'string',num2str(number_of_normal_beats))
set(handles.nOfPVC,'string',num2str(number_of_premature_ventricular_contractions))
set(handles.nOfVF,'string',num2str(number_of_ventricular_flutter_beats))
set(handles.nO2DHB,'string',num2str(number_of_second_degree_heart_block_beats))

% Additional features
handles.hrInBpm = 60./(handles.rrIntervalsInMs./1000);

% Heart rate variability (HRV) features
handles.RMSSD = sqrt( mean( (handles.rrIntervalsInMs(2:end)-handles.rrIntervalsInMs(1:end-1)).^2 ) );
handles.SDRR = std(handles.rrIntervalsInMs);
handles.SDSD = std( handles.rrIntervalsInMs(1:end-1)-handles.rrIntervalsInMs(2:end) );
handles.SD1 = (1/sqrt(2))*handles.SDSD;
handles.SD2 = sqrt( 2*handles.SDRR^2 - 0.5*handles.SDSD^2 );

RMSSD = handles.RMSSD %#ok<*NOPRT>
SDRR = handles.SDRR
SDSD = handles.SDSD
SD1 = handles.SD1
SD2 = handles.SD2


figure(2)
subplot(3,1,1)
plot(handles.rPeakTimeStamps(2:end),handles.rrIntervalsInMs,'r')
grid on
leg = legend('RR-Intervals');
set(leg,'fontsize',16)
ylabel('RRI [ms]','fontsize',16)

subplot(3,1,2)
plot(handles.rPeakTimeStamps(2:end),handles.hrInBpm ,'b')
grid on
leg = legend('HR');
set(leg,'fontsize',16)
ylabel('HR [bpm]','fontsize',16)

subplot(3,1,3)
plot(handles.rPeakTimeStamps(2:end), handles.rrCategory)
axis([0,length(handles.rrCategory),0,5])
grid on
leg = legend('Arrhythmia category');
set(leg,'fontsize',16)
ylabel('Category','fontsize',16)

figure(3)
for pp=1:length(handles.rrIntervalsInMs)-1
    if handles.rrCategory(pp) == 1
        plot(handles.rrIntervalsInMs(pp),handles.rrIntervalsInMs(pp+1),'k.')
    elseif handles.rrCategory(pp) == 2
        plot(handles.rrIntervalsInMs(pp),handles.rrIntervalsInMs(pp+1),'r.')
    elseif handles.rrCategory(pp) == 3
        plot(handles.rrIntervalsInMs(pp),handles.rrIntervalsInMs(pp+1),'m.')
    elseif handles.rrCategory(pp) == 4
        plot(handles.rrIntervalsInMs(pp),handles.rrIntervalsInMs(pp+1),'b.')
    end
    if pp == 1
        hold on
    end
end
hold off
grid on
title(sprintf('Poincare plot. RMSSD: %.1f, SDRR: %.1f, SDSD: %.1f, SD1: %.1f, SD2: %.1f.',handles.RMSSD,handles.SDRR,handles.SDSD,handles.SD1,handles.SD2))

% Update handles structure
guidata(hObject, handles);
