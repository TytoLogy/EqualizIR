function varargout = EqualizIR(varargin)
% EQUALIZIR MATLAB code for EqualizIR.fig
%      EQUALIZIR, by itself, creates a new EQUALIZIR or raises the existing
%      singleton*.
%
%      H = EQUALIZIR returns the handle to a new EQUALIZIR or the handle to
%      the existing singleton*.
%
%      EQUALIZIR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EQUALIZIR.M with the given input arguments.
%
%      EQUALIZIR('Property','Value',...) creates a new EQUALIZIR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before EqualizIR_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to EqualizIR_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Last Modified by GUIDE v2.5 04-Jan-2017 15:40:34

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @EqualizIR_OpeningFcn, ...
                   'gui_OutputFcn',  @EqualizIR_OutputFcn, ...
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

%-------------------------------------------------------------------------
% --- Executes just before EqualizIR is made visible.
%-------------------------------------------------------------------------
function EqualizIR_OpeningFcn(hObject, eventdata, handles, varargin) %#ok<*INUSL>
	%-------------------------------------------------------------
	% Initialization
	%-------------------------------------------------------------
	% Choose default command line output for EqualizIR
	handles.output = hObject;
	%-------------------------------------------------------------
	% Create internal data and defaults in handles struct
	%-------------------------------------------------------------
	handles.calfile = '';
	handles.eqfile = '';
	handles.EQ = struct(	'Fs', 100000, ...
								'NFFT', 1024, ...
								'NZ', 10, ...
								'NP', 20, ...
								'InterpMethod', 'spline', ...
								'EQMethod', 'compress', ...
								'caldata', [], ...
								'f', [], ...
								'G', [], ...
								'A', 'B', ...
								'Finv', [] ...
							);
	guidata(hObject, handles);
	%-------------------------------------------------------------
	% Update GUI from settings
	%-------------------------------------------------------------
	update_ui_str(handles.Fs_edit, handles.EQ.Fs);
	update_ui_str(handles.NFFT_edit, handles.EQ.NFFT);
	update_ui_str(handles.NZ_edit, handles.EQ.NZ);
	update_ui_str(handles.NP_edit, handles.EQ.NP);
	update_ui_val(handles.InterpMethod_popup, 1);
	update_ui_val(handles.EQMethod_popup, 1);
	%-------------------------------------------------------------
	% Update handles structure
	%-------------------------------------------------------------
	guidata(hObject, handles);
%-------------------------------------------------------------------------

%-------------------------------------------------------------------------
% --- Outputs from this function are returned to the command line.
%-------------------------------------------------------------------------
function varargout = EqualizIR_OutputFcn(hObject, eventdata, handles) 
% Get default command line output from handles structure
varargout{1} = handles.output;
%-------------------------------------------------------------------------

%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
% Load Calibration Data
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
function LoadCal_button_Callback(hObject, eventdata, handles) %#ok<*DEFNU>
	[fname, fpath] = uigetfile( {'*.cal'; '*_cal.mat'}, ...
									'Load calibration data from file...');
	if fname ~=0
		handles.calfile = fullfile(fpath, fname);	
		handles.EQ.caldata = load_headphone_cal(handles.calfile);
		plot(	handles.Cal_axes, ...
				0.001*handles.EQ.caldata.freq, ...
				handles.EQ.caldata.mag(1, :), '.-');
		ylim(handles.Cal_axes, ...
				[0.9*min(handles.EQ.caldata.mag(1, :)) ...
					1.1*max(handles.EQ.caldata.mag(1, :))]);
		grid(	handles.Cal_axes, 'on');
	end
	guidata(hObject, handles);
% 	SmoothCalCtrl_Callback(hObject, eventdata, handles);
%-------------------------------------------------------------------------

%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
% EQ Save/Load/Build
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
function SaveEQ_button_Callback(hObject, eventdata, handles)
	[fname, fpath] = uiputfile('*.eq', 'Save equalization data to file');
	if fname ~= 0
		EQ = handles.EQ; %#ok<NASGU>
		save(fullfile(fpath, fname), 'EQ', '-MAT');
	end
%-------------------------------------------------------------------------
function LoadEQ_button_Callback(hObject, eventdata, handles)
	[fname, fpath] = uigetfile( '*.eq', ...
									'Load equalization data from file...');
	if fname ~=0
		handles.eqfile = fullfile(fpath, fname);
		handles.EQ = load(handles.eqfile, '-MAT');
		plot(	handles.Cal_axes, ...
				0.001*handles.EQ.caldata.freq, ...
				handles.EQ.caldata.mag(1, :), '.-');
		ylim(handles.Cal_axes, ...
				[0.9*min(handles.EQ.caldata.mag(1, :)) ...
					1.1*max(handles.EQ.caldata.mag(1, :))]);
		grid(	handles.Cal_axes, 'on');
	end
	guidata(hObject, handles);
% 	SmoothCalCtrl_Callback(hObject, eventdata, handles);
%-------------------------------------------------------------------------

%-------------------------------------------------------------------------
function BuildEQ_button_Callback(hObject, eventdata, handles)
	%-----------------------------------------------
	% check that calibration data have been loaded
	%-----------------------------------------------
	if isempty(handles.EQ.caldata)
		error('%s: no caldata loaded!', mfilename)
	end
	%-----------------------------------------------
	% get some parameters
	%-----------------------------------------------	
	% number of gain measurements (from calibration data)
	NG = length(handles.EQ.caldata.freq);
	% frequencies from calibration data
	calfreqs = handles.EQ.caldata.freq;
	% range of frequencies
	fmin = calfreqs(1);
	fmax = calfreqs(end);
	%-----------------------------------------------	
	% compute correction xfer function
	%-----------------------------------------------	
	% process caldata gain values to convert into correction factors
	rawmags = handles.EQ.caldata.mag(1, :);
	% smooth correction
	smoothmags = sgolayfilt(raw_corr, 1, 9);
	
	% generate EQ curve according to selected method
	mtype = upper(handles.EQ.EQMethod);
	switch mtype
		case 'BOOST'
			% normalize by finding deviation from peak
			peakmag = max(calmag);
			Magnorm = peakmag - raw;		case 'ATTEN'
			COMPMETHOD = 'ATTEN';
		case 'COMPRESS'
			% shift to compromise between boost and cut
			maxdiff = max(smoothmags) - smoothmags;

			midcorr = 0.5*(max(maxdiff) - min(maxdiff));
			adj_corr = sm_corr - midcorr;
		otherwise
			fprintf('%s: unknown compensation method %s\n', ...
																	mfilename, mtype);
			fprintf('\tUsing default, BOOST method\n');
			COMPMETHOD = 'BOOST';
	end	
	
	
	if strcmpi(handles.EQ.EQMethod, 'compress')

	elseif	
	% plot xfer function
	figure(10)
	fp = f * 0.001;
	plot(fp, raw_corr, 'k', fp, sm_corr, 'r.-', fp, adj_corr, 'b.-')
	legend({'raw correction', 'smoothed correction', 'balanced correction'})
	xlabel('Frequency (kHz)');
	ylabel('Gain (db)')
	grid('on')

	
	

%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
% EQ settings
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
function Fs_edit_Callback(hObject, eventdata, handles)
	handles.EQ.Fs = read_ui_str(hObject, 'n');
	guidata(hObject, handles);
%-------------------------------------------------------------------------
function NFFT_edit_Callback(hObject, eventdata, handles)
	handles.EQ.NFFT = read_ui_str(hObject, 'n');
	guidata(hObject, handles);
%-------------------------------------------------------------------------
function NZ_edit_Callback(hObject, eventdata, handles)
	handles.EQ.NZ = read_ui_str(hObject, 'n');
	guidata(hObject, handles);
%-------------------------------------------------------------------------
function NP_edit_Callback(hObject, eventdata, handles)
	handles.EQ.NP = read_ui_str(hObject, 'n');
	guidata(hObject, handles);
%-------------------------------------------------------------------------
function InterpMethod_popup_Callback(hObject, eventdata, handles)
	% get list of strings
	mlist = cellstr(read_ui_str(hObject));
	% get selected string
	handles.EQ.InterpMethod = mlist{read_ui_val(hObject)};
	guidata(hObject, handles);
%-------------------------------------------------------------------------
function EQMethod_popup_Callback(hObject, eventdata, handles)
	% get list of strings
	mlist = cellstr(read_ui_str(hObject));
	% get selected string
	handles.EQ.EQMethod = mlist{read_ui_val(hObject)};
	guidata(hObject, handles);
%-------------------------------------------------------------------------


%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
% --- Executes during object creation, after setting all properties.
%-------------------------------------------------------------------------
function Fs_edit_CreateFcn(hObject, eventdata, handles) %#ok<*INUSD>
	if ispc && isequal(get(hObject,'BackgroundColor'), ...
								get(0,'defaultUicontrolBackgroundColor'))
		 set(hObject,'BackgroundColor','white');
	end
function NFFT_edit_CreateFcn(hObject, eventdata, handles)
	if ispc && isequal(get(hObject,'BackgroundColor'), ...
								get(0,'defaultUicontrolBackgroundColor'))
		 set(hObject,'BackgroundColor','white');
	end
function NZ_edit_CreateFcn(hObject, eventdata, handles)
	if ispc && isequal(get(hObject,'BackgroundColor'), ...
								get(0,'defaultUicontrolBackgroundColor'))
		 set(hObject,'BackgroundColor','white');
	end
function NP_edit_CreateFcn(hObject, eventdata, handles)
	if ispc && isequal(get(hObject,'BackgroundColor'), ...
								get(0,'defaultUicontrolBackgroundColor'))
		 set(hObject,'BackgroundColor','white');
	end
function InterpMethod_popup_CreateFcn(hObject, eventdata, handles)
	if ispc && isequal(get(hObject,'BackgroundColor'), ...
								get(0,'defaultUicontrolBackgroundColor'))
		 set(hObject,'BackgroundColor','white');
	end
function EQMethod_popup_CreateFcn(hObject, eventdata, handles)
	if ispc && isequal(get(hObject,'BackgroundColor'), ...
								get(0,'defaultUicontrolBackgroundColor'))
		 set(hObject,'BackgroundColor','white');
	end
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------


