
testpath = '/Users/sshanbhag/Work/Code/Matlab/dev/Audio/TimeDomainIdeas/test/';
eqfile = fullfile(testpath, 'test300.eq');
% wavfile = fullfile(testpath, 'High1_DwaveOutput_5V.WAV');
% adjfile = fullfile(testpath, 'High1_DwaveOutput_5V_eq.WAV'];
wavfile = fullfile(testpath, 'High1_CalRig_5V_unadjusted.WAV');
adjfile = fullfile(testpath, 'High1_CalRig_5V_eq.WAV');
rampms = 5;
plotspect = 'n';

ApplyFilterToWAV(	'EQFILE', eqfile, ...
						'WAVFILE', wavfile, ...
						'ADJFILE', adjfile, ...
						'RAMPMS', 6, ...
						'PLOTSPECTRA', 'n');

%%				

winfo = audioinfo(wavfile);
[wdata, Fs] = audioread(wavfile, 'native');
adjinfo = audioinfo(adjfile);
[adjdata, Fsf] = audioread(adjfile, 'native');
wdata = double(wdata);
adjdata = double(adjdata);

%%
specg_win = 1024;
specg_olap = 128;
specg_fftwin = 2048;

pbins = (2*Fs):(4*Fs);

figure(1)
subplot(211)
spectrogram(wdata(pbins), specg_win, specg_olap, specg_fftwin, Fs, 'yaxis');
colormap('gray');
map=colormap;
colormap(1-map);
colorbar
title(wavfile, 'Interpreter', 'none')

subplot(212)
spectrogram(adjdata(pbins), specg_win, specg_olap, specg_fftwin, Fsf, 'yaxis');
colormap('gray');
map=colormap;
colormap(1-map);
colorbar
title(adjfile, 'Interpreter', 'none')

%{
 
Things to note:
	(1) need to be careful that silent sections are truly 0 or else noise
	gets amplified...

%}
%%
