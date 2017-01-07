function ApplyFilterToWAV(varargin)
%------------------------------------------------------------------------
% ApplyFilterToWAV
%------------------------------------------------------------------------
% 
% given filter information in EQ struct, loads a wav file, applies filter
% and then saves to new wav file.
% 
% designed to be called from EqualizIR application!
%------------------------------------------------------------------------
% Input Arguments:
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
% 			f			[DC:Fs/2] frequencies for full equalization spectrum
% 			G			[DC:Fs/2] Gain values for equalization
% 			A, B		denominator, numerator filter coefficients that can be 
% 						used with the filter or filtfilt MATLAB commands to apply
% 						the equalization to sound arrays
% 			Finv		invert struct from josInv output
%
% Output Arguments:
% 	NONE
%------------------------------------------------------------------------
% See also: EqualizIR, filtfilt
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% Sharad J. Shanbhag
% sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: 6 January, 2017 (SJS)
%
% Revisions:
%------------------------------------------------------------------------
% TO DO:
%------------------------------------------------------------------------

%------------------------------------------
% plotting options for spectrograms
%------------------------------------------
specg_win = 512;
specg_olap = 500;
specg_fftwin = 2048;
%------------------------------------------
% If no input arg, load EQ data
%------------------------------------------
if isempty(varargin)
	[eqname, eqpath] = uigetfile( '*.eq', 'Load .EQ file...');
	if eqname == 0
		return
	else
		EQ = load(fullfile(eqpath, eqname), '-MAT', 'EQ');
	end
else
	EQ = varargin{1};
end
%------------------------------------------
% load sound (wav) file
%------------------------------------------
[callname, callpath] = uigetfile( '*.wav', 'Load .WAV file...');
if callname == 0
	return
else
	% create output name here
	[~, fname, fext] = fileparts(callname);
	adjname = [fname '_eq' fext];
end
% get file info
winfo = audioinfo(fullfile(callpath, callname));
% make sure sample rates match!!!!
if winfo.SampleRate ~= EQ.Fs
	errordlg(sprintf( ['WAV sample rate (%d) ' ...
							'does not match EQ sample rate (%d)'], ...
							winfo.SampleRate, EQ.Fs));
	return
end
% load file
fprintf('Loading file %s ... \n', fullfile(callpath, callname));
[wdata, Fs] = audioread(fullfile(callpath, callname), 'native');
fprintf('...done.\n');
% need wdata to be a row vector for some future functions to work
if ~isrow(wdata)
	wdata = wdata';
end
%------------------------------------------
% ask about onset/offset ramp duration
%------------------------------------------
rampdur = uiaskvalue('value', 5, ...
							'questiontext', 'Enter on/off ramp duration', ...
							'valuetext', 'Ramp Duration (ms)');
if rampdur < 1
	errordlg('Ramp duration must be greater than 0!');
	return
end
%------------------------------------------
% convert to double, normalize, apply ramp
%------------------------------------------
fprintf('Converting to double precision, normalizing, applying ramp\n');
wdata = sin2array(normalize(double(wdata)), rampdur, Fs);
%------------------------------------------
% then apply filter
%------------------------------------------
fprintf('Applying filter using filtfilt function and normalizing to 0.95...\n');
tic
wdataf = 0.95*normalize(filtfilt(EQ.B, EQ.A, wdata));
Tfilt = toc;
fprintf('...done in %.4f seconds\n', Tfilt);
%------------------------------------------
% plots
%------------------------------------------
%{
% plot spectra
fprintf('Plotting spectra...');
fftdbplot(wdata, Fs, 7);
fftdbplot(wdataf, Fs, 8);
fprintf('... done');
% plot spectrograms
figure(9)
% raw
subplot(211)
fprintf('Plotting raw spectrogram ...\n');
tic
spectrogram(wdata, specg_win, specg_olap, specg_fftwin, Fs, 'yaxis');
colorbar
title(fprintf('%s: raw', callname));
drawnow
Trawspec = toc;
fprintf('...done in %.4f seconds\n', Trawspec);
% filtered
subplot(212)
fprintf('Plotting adj spectrogram ...\n');
tic
spectrogram(wdataf, specg_win, specg_olap, specg_fftwin, Fs, 'yaxis');
colorbar
title(fprintf('%s: adj', adjname));
drawnow
Tadjspec = toc;
fprintf('...done in %.4f seconds\n', Tadjspec);
%}
%------------------------------------------
% Save file
%------------------------------------------
% get filename and path
[adjname, adjpath] = uiputfile(fullfile(callpath, adjname), ...
									'Save equalized WAV file');
if fname ~= 0
	% if user didn't hit cancel button, save WAV
	audiowrite(fullfile(adjpath, adjname), wdataf, Fs, ...
						'BitsPerSample', winfo.BitsPerSample);
end


