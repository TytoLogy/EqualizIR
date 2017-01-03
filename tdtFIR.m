function [filtcoefs, varargout] = tdtFIR(freqs, norms, varargin)
%------------------------------------------------------------------------
% [filtcoefs, varargout] = tdtFIR(freqs, norms, ntaps)
%------------------------------------------------------------------------
% 
% computes FIR filter from input freqs and normalization
% 
%------------------------------------------------------------------------
% Input Arguments:
%
%	freqs		frequencies (in range [0 1] where 1 corresponds to Fs
%	norms		normalization coefficients in linear (non-dB) scale
%
%	Optional:
%		ntaps			# of taps for output filter. default is 1024
% 
% Output Arguments:
%	filtcoefs		[ntaps X 1] array of coefficients
%------------------------------------------------------------------------
% See also: fir2
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% Sharad J. Shanbhag
% sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: 15 November, 2016 (SJS)
%	based on RZ6SigCalFIR.m from TDT
%
% Revisions:
%------------------------------------------------------------------------
% TO DO:
%------------------------------------------------------------------------

ntaps = 1024;

if nargin == 3
	ntaps = varargin{1};
end

% make sure beginning is zero
if freqs(1) ~= 0
	freqs = [0 freqs];
	norms = [norms(1) norms];
end

% make sure end is 1;
if freqs(end) ~= 1
	freqs = [freqs 1];
	norms = [norms norms(end)];
end

% calculate filter coefficients
filtcoefs = fir2(ntaps, freqs, norms);

if nargout > 1
	varargout{1}.freqs = freqs;
	varargout{1}.norms = norms;
end
