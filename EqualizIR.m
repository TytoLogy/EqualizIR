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

% Last Modified by GUIDE v2.5 09-Jan-2017 17:44:43

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
	%
	% Program-defined fields of handles struct:
	% 		calfile	name of calibration data file
	% 		eqfile	name of equalization data file
	% 		EQ		EQualization data struct:
	% 			Fs			sample rate for impulse response, filter (samples/sec)
	% 							¡¡this must match that of sound files to be processed!!!
	% 			NFFT		length of impulse response (samples)
	% 			NZ			number of zeros for filter
	% 			NP			number of poles for filter
	% 			InterpMethod	Interpolation method used to complete correction
	% 								spectrum (sent as argument to josInv.m function)
	% 			EQMethod			method used to compute equalization factor
	% 									{'COMPRESS', 'ATTEN', 'BOOST'}
	% 			TargetLevel		level used as target by COMPRESS method
	% 			caldata	calibration data struct
	%			CalSmoothMethod	method used to smooth calibration data before
	% 									generating filter
	% 										{'SAVGOL', 'MOVAVG'}
	% 			CalSmoothParameters	settings for smoothing method
	% 										For SAVGOL, [order, framesize]
	% 										For MOVAVG, [window size]
	% 			EQmags	equalization values calculated by EQ method
	% 							¡calculated at frequencies specified in caldata struct, 
	% 							 NOT at frequencies in f, G arrays - those are calculated
	% 							 separately using the EQmags!
	%			CorrectionLimit	limit to correction value in dB. if 0, no limit!
	%			FreqLimit	limit to frequency range 
	%								(0 if no limit, 1 if limit enabled)
	%			FLimitMin		initially set to caldata.freq(1)
	%			FLimitMax		initially set to caldata.freq(end)
	% 			f			[DC:Fs/2] frequencies for full equalization spectrum
	% 			G			[DC:Fs/2] Gain values for equalization
	% 			A, B		denominator, numerator filter coefficients that can be 
	% 						used with the filter or filtfilt MATLAB commands to apply
	% 						the equalization to sound arrays
	% 			Finv		invert struct from josInv output
	%-------------------------------------------------------------
	handles.calfile = '';
	handles.eqfile = '';
	handles.EQ = struct(	'Fs', 250000, ...
								'NFFT', 1024, ...
								'NZ', 10, ...
								'NP', 20, ...
								'InterpMethod', 'linear', ...
								'EQMethod', 'compress', ...
								'TargetLevel', 80, ...
								'caldata', [], ...
								'CalSmoothMethod', 'SAVGOL', ...
								'CalSmoothParameters', [], ...
								'EQmags', [], ...
								'CorrectionLimit', 0, ...
								'FreqLimit', 0, ...
								'FLimitMin', 0, ...
								'FLimitMax', 125000, ...
								'R', [], ...
								'f', [], ...
								'G', [], ...
								'A', [], ...
								'B', [], ...
								'Finv', [] ...
							);
	% need to assign this outside of struct() in order to prevent Matlab
	% from creating EQ and a 1X2 struct array!
	handles.EQ.CalSmoothParameters = {[1 9], 5};
	guidata(hObject, handles);
	%-------------------------------------------------------------
	% Update GUI from settings
	%-------------------------------------------------------------
	update_ui_str(handles.Fs_edit, handles.EQ.Fs);
	update_ui_str(handles.NFFT_edit, handles.EQ.NFFT);
	update_ui_str(handles.NZ_edit, handles.EQ.NZ);
	update_ui_str(handles.NP_edit, handles.EQ.NP);
	update_ui_val(handles.InterpMethod_popup, 2);
	update_ui_val(handles.EQMethod_popup, 1);
	% set MiddleLevel_radiobutton value to 1 (selected)
	update_ui_val(handles.MiddleLevel_radiobutton, 1);
	% set TargetLevel_radiobutton value to 0 (unselected)
	update_ui_val(handles.TargetLevel_radiobutton, 0);
	% update target level and middle level text boxes
	update_ui_str(handles.TargetLevel_edit, handles.EQ.TargetLevel);
	update_ui_str(handles.MiddleLevel_text, '--');
	% since default is compress, middle level, disable the target level edit
	% box
	disable_ui(handles.TargetLevel_edit);
	% update correction limit
	if handles.EQ.CorrectionLimit
		update_ui_val(handles.CorrectionLimit_checkbox, 1);
		update_ui_str(handles.CorrectionLimit_edit, handles.EQ.CorrectionLimit);
		enable_ui([handles.CorrectionLimit_text, handles.CorrectionLimit_edit]);
	else
		update_ui_val(handles.CorrectionLimit_checkbox, 0);
		disable_ui([handles.CorrectionLimit_text, handles.CorrectionLimit_edit]);
	end
	% freq range
	update_ui_val(handles.LimitFreqRange_checkbox, handles.EQ.FreqLimit);
	update_ui_str(handles.FLimitMin_edit, handles.EQ.FLimitMin);
	update_ui_str(handles.FLimitMax_edit, handles.EQ.FLimitMax);
	% smoothing settings
	switch upper(handles.EQ.CalSmoothMethod)
		case 'SAVGOL'
			enable_ui(	[handles.SmoothVal1_text, handles.SmoothVal1_edit, ...
								handles.SmoothVal2_text, handles.SmoothVal2_edit]);
			update_ui_str(handles.SmoothVal1_text, 'order');
			update_ui_str(handles.SmoothVal1_edit, ...
									handles.EQ.CalSmoothParameters{1}(1));
			update_ui_str(handles.SmoothVal2_text, 'frame');
			update_ui_str(handles.SmoothVal2_edit, ...
									handles.EQ.CalSmoothParameters{1}(2));
		case 'MOVAVG'
			enable_ui([handles.SmoothVal1_text, handles.SmoothVal1_edit]);
			disable_ui([handles.SmoothVal2_text, handles.SmoothVal2_edit]);
			update_ui_str(handles.SmoothVal1_text, 'window');
			update_ui_str(handles.SmoothVal1_edit, ...
									handles.EQ.CalSmoothParameters{2});
	end
	% clear plots
	cla(handles.Cal_axes);
	cla(handles.EQ_axes);
	cla(handles.Filter_axes);
	cla(handles.Desired_axes);
	cla(handles.IR_axes);
	cla(handles.PZ_axes);
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
% Calibration Data Functions
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
function SmoothCal_ctrl_Callback(hObject, eventdata, handles) %#ok<*DEFNU>
	% get value of smooth method
	smoothmethod = read_ui_val(hObject);
	switch(smoothmethod)
		case 1
			handles.EQ.CalSmoothMethod = 'SAVGOL';
			enable_ui(	[handles.SmoothVal1_text, handles.SmoothVal1_edit, ...
								handles.SmoothVal2_text, handles.SmoothVal2_edit]);
			update_ui_str(handles.SmoothVal1_text, 'order');
			update_ui_str(handles.SmoothVal1_edit, ...
									handles.EQ.CalSmoothParameters{1}(1));
			update_ui_str(handles.SmoothVal2_text, 'frame');
			update_ui_str(handles.SmoothVal2_edit, ...
									handles.EQ.CalSmoothParameters{1}(2));
		case 2
			handles.EQ.CalSmoothMethod = 'MOVAVG';
			enable_ui([handles.SmoothVal1_text, handles.SmoothVal1_edit]);
			disable_ui([handles.SmoothVal2_text, handles.SmoothVal2_edit]);
			update_ui_str(handles.SmoothVal1_text, 'window');
			update_ui_str(handles.SmoothVal1_edit, ...
									handles.EQ.CalSmoothParameters{2});
	end
	guidata(hObject, handles);
%-------------------------------------------------------------------------
function SmoothVal1_edit_Callback(hObject, eventdata, handles)
	tmp = read_ui_str(handles.SmoothVal1_edit, 'n');
	method = read_ui_val(handles.SmoothCal_ctrl);
	if tmp > 0
		handles.EQ.CalSmoothParameters{method}(1) = tmp;
	else
		update_ui_str(handles.SmoothVal1_edit, ...
									handles.EQ.CalSmoothParameters{method}(1));
		errordlg('Value must be greater than 0', 'EqualizIR Error');
	end
	guidata(hObject, handles)
%-------------------------------------------------------------------------
function SmoothVal2_edit_Callback(hObject, eventdata, handles)
	tmp = read_ui_str(handles.SmoothVal2_edit, 'n');
	if tmp > 0
		if even(tmp)
			update_ui_str(handles.SmoothVal2_edit, ...
								handles.EQ.CalSmoothParameters{1}(2));
			errordlg('Value must be odd', 'EqualizIR Error');
			return
		else
			handles.EQ.CalSmoothParameters{1}(2) = tmp;
		end
	else
		update_ui_str(handles.SmoothVal2_edit, ...
									handles.EQ.CalSmoothParameters{1}(2));
		errordlg('Value must be greater than 0', 'EqualizIR Error');
	end
	guidata(hObject, handles)
%-------------------------------------------------------------------------
function LimitFreqRange_checkbox_Callback(hObject, eventdata, handles)
	%-----------------------------------------------
	% checks
	%-----------------------------------------------
	%  see if calibration data have been loaded
	if isempty(handles.EQ.caldata)
		errordlg('%no caldata loaded!');
		return
	end
	%-----------------------------------------------
	% update things
	%-----------------------------------------------
	handles.EQ.FreqLimit = read_ui_val(hObject);
	if handles.EQ.FreqLimit
		handles.EQ.FLimitMin = read_ui_str(handles.FLimitMin_edit, 'n');
		handles.EQ.FLimitMax = read_ui_str(handles.FLimitMax_edit, 'n');
	end
	guidata(hObject, handles);
%---------------------------------------------------------------------
function FLimitMin_edit_Callback(hObject, eventdata, handles)
	%-----------------------------------------------
	% checks
	%-----------------------------------------------
	%  see if calibration data have been loaded
	if isempty(handles.EQ.caldata)
		errordlg('%no caldata loaded!');
		return
	end
	%-----------------------------------------------
	% update things
	%-----------------------------------------------
	handles.EQ.FLimitMin = read_ui_str(handles.FLimitMin_edit, 'n');
	guidata(hObject, handles);
%---------------------------------------------------------------------
function FLimitMax_edit_Callback(hObject, eventdata, handles)
	%-----------------------------------------------
	% checks
	%-----------------------------------------------
	%  see if calibration data have been loaded
	if isempty(handles.EQ.caldata)
		errordlg('%no caldata loaded!');
		return
	end
	%-----------------------------------------------
	% update things
	%-----------------------------------------------
	handles.EQ.FLimitMax = read_ui_str(handles.FLimitMax_edit, 'n');
	guidata(hObject, handles);
%---------------------------------------------------------------------
function BuildEQ_button_Callback(hObject, eventdata, handles)
%---------------------------------------------------------------------
%---------------------------------------------------------------------
% Builds the equalization curve from calibration data
%	See compensate_signal.m in TytoLogy AudioToolbox for original
%	implemenation of this code
%---------------------------------------------------------------------
	%-----------------------------------------------
	% checks
	%-----------------------------------------------
	%  see if calibration data have been loaded
	if isempty(handles.EQ.caldata)
		errordlg('No calibration data are loaded!', mfilename)
		return
	end
	%-----------------------------------------------
	% get some parameters
	%-----------------------------------------------	
	% frequencies from calibration data
	calfreqs = handles.EQ.caldata.freq;
	%---------------------------------------------------------------------
	% pre-process caldata gain values to convert into correction factors
	% -> smooth correction using golay filter or moving average
	%---------------------------------------------------------------------
	if strcmpi(handles.EQ.CalSmoothMethod, 'SAVGOL')
		smoothmags = sgolayfilt(handles.EQ.caldata.mag(1, :), ...
										handles.EQ.CalSmoothParameters{1}(1), ...
										handles.EQ.CalSmoothParameters{1}(2));
	else
		smoothmags = moving_average(handles.EQ.caldata.mag(1, :), ...
										handles.EQ.CalSmoothParameters{2});
	end
	% store smoothed values
	handles.smoothmags = smoothmags;
	guidata(hObject, handles);
	% plot smoothed values
	hold(handles.Cal_axes, 'on');
		plot(handles.Cal_axes, 0.001*calfreqs, smoothmags, 'r-');
	hold(handles.Cal_axes, 'off');
	legend(handles.Cal_axes, {'Raw', 'Smoothed'});
	%---------------------------------------------------------------------
	% limit range if necessary
	%---------------------------------------------------------------------
	if handles.EQ.FreqLimit
		% need to find max, min of calibration range
		R = find(between(	calfreqs, ...
											handles.EQ.FLimitMin, ...
											handles.EQ.FLimitMax) ...
								==1);
	else
		% use full calibration range
		R = 1:length(calfreqs);
	end
	%---------------------------------------------------------------------
	% generate EQ curve according to selected method
	%---------------------------------------------------------------------
	switch upper(handles.EQ.EQMethod)
		%------------------------------------------------------------------------
		% apply correction using BOOST method
		%------------------------------------------------------------------------
		% procedure:	find additive compensation values for frequency range 
		%					for which there are calibration data and apply to FFT, 
		%					then iFFT to get corrected version
		%------------------------------------------------------------------------
		case 'BOOST'
			% find peak magnitude, then calculate equalization values
			% by finding deviation from peak
			peakmag = max(smoothmags(R));
			handles.EQ.EQmags = peakmag - smoothmags(R);
			if handles.EQ.CorrectionLimit
				handles.EQ.EQmags = limit_correction(handles.EQ.EQmags, ...
																handles.EQ.CorrectionLimit);
			end
			guidata(hObject, handles);
		%------------------------------------------------------------------------
		% apply correction using ATTEN method
		%------------------------------------------------------------------------
		% procedure:	find subtractive compensation values for frequency range 
		%					for which there are calibration data and apply to FFT, 
		%					then iFFT to get corrected version
		%------------------------------------------------------------------------
		case 'ATTEN'
			% find lowest magnitude, then calculate equalization values
			% by finding deviation from minimum
			minmag = min(smoothmags(R));
			handles.EQ.EQmags = minmag - smoothmags(R);
			if handles.EQ.CorrectionLimit
				handles.EQ.EQmags = limit_correction(handles.EQ.EQmags, ...
																handles.EQ.CorrectionLimit);
			end
			guidata(hObject, handles);
		%------------------------------------------------------------------------
		% apply correction using COMPRESS method
		%------------------------------------------------------------------------
		% procedure:	find subtractive compensation values for frequency range 
		%					for which there are calibration data and apply to FFT, 
		%					then iFFT to get corrected version
		%------------------------------------------------------------------------
		% some assumptions:
		% 					magnitude values are in ACTUAL, dB SPL range.  
		% 						¡this algorithm blows up for negative magnitudes!
		%------------------------------------------------------------------------
		case 'COMPRESS'
			% see if Middle or Target level is being used
			if read_ui_val(handles.MiddleLevel_radiobutton)
				%------------------------------------------------------
				% Use Middle level as target
				%  This is a compromise between boost and atten 
				%  by finding middle of dB range
				%------------------------------------------------------
					% find max and min in magnitude spectrum
				maxmag = max(smoothmags(R));
				minmag = min(smoothmags(R));
				% compute middle value
				midmag = ((maxmag - minmag) / 2) + minmag;
				% normalize by finding deviation from middle level
				handles.EQ.EQmags = midmag - smoothmags(R);			
				% update GUI
				update_ui_str(handles.MiddleLevel_text, handles.EQ.TargetLevel);
			else
				%------------------------------------------------------
				% Use specified Target level to compute boost/atten
				%------------------------------------------------------
				handles.EQ.TargetLevel = read_ui_val(handles.EQ.TargetLevel_edit);
				handles.EQ.EQmags = handles.EQ.TargetLevel - smoothmags(R);
			end
			% Limit correction if specified
			if handles.EQ.CorrectionLimit
				handles.EQ.EQmags = limit_correction(handles.EQ.EQmags, ...
															handles.EQ.CorrectionLimit);
			end
			guidata(hObject, handles);
		otherwise
			errordlg(sprintf('unsupported compensation method %s\n', ...
															handles.EQ.EQMethod));
			return
	end
	% plot xfer function, storing handles to plots in handles.Corr_H
	plot(handles.EQ_axes, 0.001*calfreqs(R), handles.EQ.EQmags);
	legend(handles.EQ_axes, sprintf('%s EQ', lower(handles.EQ.EQMethod)));
	ylabel(handles.EQ_axes, 'Gain (db)')
	grid(handles.EQ_axes, 'on')
	box(handles.EQ_axes, 'off');
	xlim(handles.EQ_axes, 0.001*[min(calfreqs) max(calfreqs)]);
% 	set(handles.EQ_axes, 'XTickLabel', '');
	handles.EQ.R = R;
	guidata(hObject, handles);
%-------------------------------------------------------------------------

%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
% EQ Build
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
function BuildFilter_button_Callback(hObject, eventdata, handles)
	%-----------------------------------------------
	% checks
	%-----------------------------------------------
	%  see if calibration data have been loaded
	if isempty(handles.EQ.caldata)
		errordlg('No Calibration Data are Loaded!');
		return
	end
	% see if eq curve is built
	if isempty(handles.EQ.EQmags)
		errordlg('No Correction data yet - try pressing Build Correction!')
		return
	end	
	% make sure Fs is going to work with calibration data
	if handles.EQ.Fs/2 < handles.EQ.caldata.freq(end)
		errordlg(sprintf(['Signal Nyquist frequency (%d) is lower than ' ...
				 'max calibration frequency (%d)']), ...
				 handles.EQ.Fs/2, handles.EQ.caldata.freq(end));
	end
	%-----------------------------------------------
	% get some parameters
	%-----------------------------------------------	
	% frequencies from calibration data
	calfreqs = handles.EQ.caldata.freq(handles.EQ.R);
	%------------------------------------------------------------------------
	% some final tidying and ...
	%------------------------------------------------------------------------
	% Gain measurements - make a working copy of EQ curve
	Gdb = handles.EQ.EQmags;
	% number of gain measurements (from calibration data)
	NG = length(Gdb);
	% Must decide on a dc value.
	% Either use what is known to be true or pick something "maximally
	% smooth".  Here we do a simple linear extrapolation:
	dc_amp = Gdb(1) - calfreqs(1)*(Gdb(2)-Gdb(1))/(calfreqs(2)-calfreqs(1));
	% JOS: (original code)
	% Must also decide on a value at half the sampling rate.
	% Use either a realistic estimate or something "maximally smooth".
	% Here we do a simple linear extrapolation. While zeroing it
	% is appealing, we do not want any zeros on the unit circle here.
	Gdb_last_slope = (Gdb(NG) - Gdb(NG-1)) / (calfreqs(NG) - calfreqs(NG-1));
	% SJS:
	% one problem with this approach is that if (1) the max value of measured
	% frequencies (f(NG)) is far below fs/2, and if (2) the  value of the last
	% slope, Gdb_last_slope, is positive, that final value for nyq_amp is
	% going to "blow up" positively.
	% --- a check on this .... ----
	% try checking the last slope value and, if positive, do something else...
	if Gdb_last_slope > 0
		nadjpts = 10;
		adjindx = (NG-nadjpts+1):NG;
		x = linspace(0, 1, nadjpts);
		% use neg. parabolic to calculate values
		y = 1 - .15*x.^2;
		% apply correction factor...
		Gdb(adjindx) = y.*Gdb(adjindx);
		Gdb_last_slope = (Gdb(NG) - Gdb(NG-1)) / (calfreqs(NG) - calfreqs(NG-1));
		% ...and plot corrected gain in EQ_axes
		hold(handles.EQ_axes, 'on');
			plot(handles.EQ_axes, 0.001*calfreqs, Gdb, 'g.-')
		hold(handles.EQ_axes, 'off');
		legend(handles.EQ_axes, {'original', 'adjusted'});
		%{
		% FOR DEBUGGING
		% plot details of correction factor
		figure(2)
		plot(x, Gdb(adjindx), 'o')
		hold on
			plot(x, y, '.k');
			plot(x, y.*Gdb(adjindx), '.r')
		hold off
		grid
		legend({'original', 'adjusted'});
		%}
	end
	% now, compute amplitude at Nyquist freq (if necessary)
	if handles.EQ.Fs / 2 ~= calfreqs(NG)
		nyq_amp = Gdb(NG) + Gdb_last_slope * (handles.EQ.Fs/2 - calfreqs(NG));
		handles.EQ.G = [dc_amp, Gdb, nyq_amp];
		handles.EQ.f = [0,calfreqs,handles.EQ.Fs/2];
	else
		handles.EQ.G = [dc_amp, Gdb];
		handles.EQ.f = [0, calfreqs];
	end
	guidata(hObject, handles);
	%------------------------------------------------------------------------
	% ...compute filter
	%------------------------------------------------------------------------
	[handles.EQ.B, handles.EQ.A, handles.EQ.Finv] = ...
				josInv(	handles.EQ.f, ...
							handles.EQ.G, ...
							handles.EQ.NZ, ...
							handles.EQ.NP, ...
							handles.EQ.NFFT, ...
							handles.EQ.Fs, ...
							'ShowPlot', 'n', ...
							'InterpMethod', handles.EQ.InterpMethod);
	guidata(hObject, handles);
	% mean squared error between desired and computed
	MSE = mean((handles.EQ.Finv.Gdbfk' - db(handles.EQ.Finv.Hh)).^2);
	update_ui_str(handles.MSE_text, sprintf('%.4f', MSE));
	
	%------------------------------------------------------------------------
	% plots
	%------------------------------------------------------------------------
	% plot Filter
	plot(	handles.Filter_axes, ...
				handles.EQ.Finv.Fk, ...
				db([handles.EQ.Finv.Smpp(:), handles.EQ.Finv.Hh(:)]));
	grid(handles.Filter_axes, 'on');
	xlabel(handles.Filter_axes, 'Frequency (Hz)');
	ylabel(handles.Filter_axes, 'Correction Magnitude (dB)');
	legend(handles.Filter_axes, 'Desired','Filter');
	xlim(handles.Filter_axes, [min(handles.EQ.Finv.Fk) max(handles.EQ.Finv.Fk)]);

	% plot measured and fit magnitude response
	semilogx(handles.Desired_axes, ...
						handles.EQ.Finv.Fk(2:end-1), ...
						handles.EQ.Finv.Gdbfk(2:end-1),'.k'); 
	hold(handles.Desired_axes, 'on'); 
		semilogx(handles.Desired_axes, handles.EQ.f, handles.EQ.G, 'o');
	hold(handles.Desired_axes, 'off');
	grid(handles.Desired_axes, 'on');
	axis(handles.Desired_axes, ...
				[	min(handles.EQ.f)/2 max(handles.EQ.f)*2 ...
					min(handles.EQ.Finv.Gdbfk) 1.1*max(handles.EQ.Finv.Gdbfk)]);
	xlabel(handles.Desired_axes, 'Frequency (Hz)');
	ylabel(handles.Desired_axes, 'Magnitude (dB)');
	legend(handles.Desired_axes, {'extra/inter/resampled', 'desired'});
	
	% plot impulse response
	plot(handles.IR_axes, handles.EQ.Finv.s(1:handles.EQ.NFFT/2), '-k');
	grid(handles.IR_axes, 'on');
	xlabel(handles.IR_axes, 'Samples');
	ylabel(handles.IR_axes, 'Amplitude');
	xlim(handles.IR_axes, [0 (0.5*handles.EQ.NFFT)+1]);

	% plot poles and zeros
	cla(handles.PZ_axes);
	zplane(handles.EQ.B, handles.EQ.A, handles.PZ_axes);
%  	axis(handles.PZ_axes, 'square')
	xlim(handles.PZ_axes, [-1.1 1.1]);
	ylim(handles.PZ_axes, [-1.1 1.1]);
	legend(handles.PZ_axes, {'zeros', 'poles'})
	
%-------------------------------------------------------------------------



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
	% should be {compress, atten, boost}
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
	switch upper(handles.EQ.EQMethod)
		case 'COMPRESS'
			enable_ui([handles.MiddleLevel_radiobutton, ...
							handles.TargetLevel_radiobutton]);
			if read_ui_val(handles.MiddleLevel_radiobutton)
				disable_ui(handles.TargetLevel_edit);
			else
				enable_ui(handles.TargetLevel_edit);
			end
		case {'ATTEN', 'BOOST'}
			disable_ui([handles.MiddleLevel_radiobutton, ...
							handles.TargetLevel_radiobutton]);
	end	
%-------------------------------------------------------------------------
function MiddleLevel_radiobutton_Callback(hObject, eventdata, handles)
% enable when "COMPRESS" eq method is selected in EQmethod_popup
% when Value is 1 (selected), disable TargetLevel_edit 
	val = read_ui_val(handles.MiddleLevel_radiobutton);
	fprintf('MiddleLevel_radiobutton clicked, value = %d\n', val);
	if val
		update_ui_val(handles.TargetLevel_radiobutton, 0);
		disable_ui(handles.TargetLevel_edit);
% 		update_ui_str(handles.MiddleLevel_text, handles.EQ.TargetLevel);
	else
		update_ui_val(handles.TargetLevel_radiobutton, 1);
		enable_ui(handles.TargetLevel_edit);
% 		update_ui_str(handles.TargetLevel_edit, handles.EQ.TargetLevel);
	end
	guidata(hObject, handles);
%-------------------------------------------------------------------------
function TargetLevel_radiobutton_Callback(hObject, eventdata, handles)
% enable, make visible when "COMPRESS" eq method is selected in
% EQmethod_popup
% when Value is 1 (selected), enable TargetLevel_edit 
	val = read_ui_val(handles.TargetLevel_radiobutton);
	fprintf('TargetLevel_radiobutton clicked, value = %d\n', val);
	if val
		update_ui_val(handles.MiddleLevel_radiobutton, 0);
		enable_ui(handles.TargetLevel_edit);
		handles.EQ.TargetLevel = round(read_ui_str(handles.TargetLevel_edit, 'n'));
	else
		update_ui_val(handles.MiddleLevel_radiobutton, 1);
		disable_ui(handles.TargetLevel_edit);
		update_ui_str(handles.MiddleLevel_text, '');
	end
	guidata(hObject, handles);
%-------------------------------------------------------------------------
function TargetLevel_edit_Callback(hObject, eventdata, handles)
% set EQ target level
	handles.EQ.TargetLevel = round(read_ui_str(handles.TargetLevel_edit, 'n'));
	guidata(hObject, handles);
%-------------------------------------------------------------------------
function CorrectionLimit_checkbox_Callback(hObject, eventdata, handles)
% if checked, enable CorrectionLimit_text, and CorrectionLimit_edit, update
% value
	checkval = read_ui_val(hObject);
	if checkval
		enable_ui(handles.CorrectionLimit_edit);
		enable_ui(handles.CorrectionLimit_text);
		handles.EQ.CorrectionLimit = ...
							read_ui_str(handles.CorrectionLimit_edit, 'n');
	else
		disable_ui(handles.CorrectionLimit_edit);
		disable_ui(handles.CorrectionLimit_text);
		handles.EQ.CorrectionLimit = 0;
	end
	guidata(hObject, handles);
%-------------------------------------------------------------------------
function CorrectionLimit_edit_Callback(hObject, eventdata, handles)
% set Correction Limit value
	handles.EQ.CorrectionLimit = read_ui_str(hObject, 'n');
	guidata(hObject, handles);
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------

%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
% UTILITY FUNCTIONS (internal)	
%-------------------------------------------------------------------------
function magsout = limit_correction(magsin, corrlimit)
% given input magnitudes, magsin, and correction limit, corrlimit), 
% limit_correction will set values of magsin outside the range 
%[-corrlimit, corrlimit] to corrlimit (-corrlimit if below range, 
% + corrlimit if above)
	if all(magsin < corrlimit)
		warning('EqualizIR:CORRLIMIT', 'CORRLIMIT > all values');
		fprintf('\tConsider raising the Target SPL level\n');
		fprintf('\tin order to balance correction amount!\n\n');
	elseif all(magsin > corrlimit)
		warning('EqualizIR:CORRLIMIT', 'CORRLIMIT < all values');
		fprintf('\tConsider lowering the Target SPL level\n');
		fprintf('\tin order to balance correction amount!\n\n');
	end
	magsout = magsin;
	magsout(magsout > corrlimit) = corrlimit;
	magsout(magsout < -corrlimit) = -corrlimit;
%-------------------------------------------------------------------------
function Debug_button_Callback(hObject, eventdata, handles)
	keyboard
%-------------------------------------------------------------------------

%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
% MENU CALLBACKS
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
function LoadCal_menu_Callback(hObject, eventdata, handles)
%---------------------------------------------
% Load, plot calibration data
%---------------------------------------------
	[fname, fpath] = uigetfile( {'*.cal'; '*_cal.mat'}, ...
									'Load calibration data from file...');
	if fname ~=0
		handles.calfile = fullfile(fpath, fname);	
		handles.EQ.caldata = load_headphone_cal(handles.calfile);
		plot(handles.Cal_axes, ...
				0.001*handles.EQ.caldata.freq, ...
				handles.EQ.caldata.mag(1, :), '.-');
		ylim(handles.Cal_axes, ...
				[0.9*min(handles.EQ.caldata.mag(1, :)) ...
					1.1*max(handles.EQ.caldata.mag(1, :))]);
		xlim(handles.Cal_axes, 0.001*[min(handles.EQ.caldata.freq) ...
												max(handles.EQ.caldata.freq)]);
		grid(handles.Cal_axes, 'on');
		ylabel(handles.Cal_axes, 'dB (SPL)')
		xlabel(handles.Cal_axes, 'Frequency (kHz)')
		box(handles.Cal_axes, 'off');
		update_ui_str(handles.FLimitMin_edit, min(handles.EQ.caldata.freq));
		update_ui_str(handles.FLimitMax_edit, max(handles.EQ.caldata.freq));
	end
	update_ui_str(handles.CalFileName_text, handles.calfile);
	guidata(hObject, handles);
%-------------------------------------------------------------------------
function EQ_menu_Callback(hObject, eventdata, handles)
%---------------------------------------------
% placeholder
%---------------------------------------------
	return
%-------------------------------------------------------------------------
function SaveEQ_menu_Callback(hObject, eventdata, handles)
%---------------------------------------------
% saves EQ data to .eq file (MAT file format)
%---------------------------------------------
	% check to make sure EQ data exist...
	if isempty(handles.EQ.G)
		% if not, throw error
		errordlg('cannot save EQ data - no data to save!');
		return
	end
	% get filename and path
	[fname, fpath] = uiputfile('*.eq', 'Save equalization data to file');
	if fname ~= 0
		% if user didn't hit cancel button, save EQ struct in mat file
		EQ = handles.EQ; %#ok<NASGU>
		save(fullfile(fpath, fname), 'EQ', '-MAT');
	end
%-------------------------------------------------------------------------
function LoadEQ_menu_Callback(hObject, eventdata, handles)
%---------------------------------------------
% loads EQ data from .eq file (MAT file format)
%---------------------------------------------
	[fname, fpath] = uigetfile( '*.eq', ...
									'Load equalization data from file...');
	if fname ~=0
		handles.eqfile = fullfile(fpath, fname);
		handles.EQ = load(handles.eqfile, '-MAT');
		plot(handles.Cal_axes, ...
				0.001*handles.EQ.caldata.freq, ...
				handles.EQ.caldata.mag(1, :), '.-');
		ylim(handles.Cal_axes, ...
				[0.9*min(handles.EQ.caldata.mag(1, :)) ...
					1.1*max(handles.EQ.caldata.mag(1, :))]);
		grid(handles.Cal_axes, 'on');
		ylabel(handles.Cal_axes, 'dB (SPL)')
		xlabel(handles.Cal_axes, 'Frequency (kHz)')
	end
	guidata(hObject, handles);
%-------------------------------------------------------------------------
%--------------------------------------------------------------------
function WAVapply_menu_Callback(hObject, eventdata, handles)
%---------------------------------------------
% applies filter to WAV file
%---------------------------------------------
	ApplyFilterToWAV('EQ', handles.EQ)
%--------------------------------------------------------------------
%--------------------------------------------------------------------
%--------------------------------------------------------------------

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
function TargetLevel_edit_CreateFcn(hObject, eventdata, handles)
	if ispc && isequal(get(hObject,'BackgroundColor'), ...
								get(0,'defaultUicontrolBackgroundColor'))
		 set(hObject,'BackgroundColor','white');
	end
function CorrectionLimit_edit_CreateFcn(hObject, eventdata, handles)
	if ispc && isequal(get(hObject,'BackgroundColor'), ...
								get(0,'defaultUicontrolBackgroundColor'))
		 set(hObject,'BackgroundColor','white');
	end
function SmoothCal_ctrl_CreateFcn(hObject, eventdata, handles)
	if ispc && isequal(get(hObject,'BackgroundColor'), ...
								get(0,'defaultUicontrolBackgroundColor'))
		 set(hObject,'BackgroundColor','white');
	end
function SmoothVal1_edit_CreateFcn(hObject, eventdata, handles)
	if ispc && isequal(get(hObject,'BackgroundColor'), ...
								get(0,'defaultUicontrolBackgroundColor'))
		 set(hObject,'BackgroundColor','white');
	end
function SmoothVal2_edit_CreateFcn(hObject, eventdata, handles)
	if ispc && isequal(get(hObject,'BackgroundColor'), ...
								get(0,'defaultUicontrolBackgroundColor'))
		 set(hObject,'BackgroundColor','white');
	end
function FLimitMin_edit_CreateFcn(hObject, eventdata, handles)
	if ispc && isequal(get(hObject,'BackgroundColor'), ...
								get(0,'defaultUicontrolBackgroundColor'))
		 set(hObject,'BackgroundColor','white');
	end
function FLimitMax_edit_CreateFcn(hObject, eventdata, handles)
	if ispc && isequal(get(hObject,'BackgroundColor'), ...
								get(0,'defaultUicontrolBackgroundColor'))
		 set(hObject,'BackgroundColor','white');
	end
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
