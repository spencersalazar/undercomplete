function [y_re] = undercomplete2(y, B, block_size, hop_size, window, weight_pole)

% function [y_re] = undercomplete2(y, B, block_size, hop_size, window, 
% 	weight_pole)
% 
% Model an input signal y using only a set of basis vectors defined by B, 
% L2+L1 optimization is used to determine the optimal weighting for the given 
% basis.
%
% This variant minimizes variation in magnitude frequency response rather 
% than in the time domain.  
% 
% === REQUIRED ARGUMENTS ===
% y
%	the input signal (column vector).
% 
% B
%	basis matrix. Each row of B is treated as the i'th basis vector. 
%
% === OPTIONAL ARGUMENTS ===
% block_size
%	length of block to operate on; should typically equal column length of B.
%	(default: 1024)
%
% hop_size
%	number of samples to advance on each modeling step
%	(default: block_size/2)
%
% window
%	block_size length column vector with the desired block window
%	(default: hanning)
% 
% weight_pole
% 	Pole location of one-pole filter applied to time-varying basis weights. 
% 	(default: 0.65)
% 
% === RETURN VALUE ===
% y_re
%	Reconstructed output signal.
%


if nargin < 2
	error('undercomplete: must provide input signal y and basis matrix B');
end

if nargin < 3
	block_size = 1024;
end
if nargin < 4
	hop_size = block_size/2;
end
if nargin < 5
	window = hanning(block_size);
end
if nargin < 6
	weight_pole = 0.65;
end

addpath('GPSR_6.0');
addpath('lasso');

weight_a1 = weight_pole;
weight_b0 = 1-weight_a1;

Jb = size(B,1);

Bfft = abs(fft(B'));

% pad out to hop size
y = [y; zeros(hop_size - mod(length(y),hop_size),1)];
y_re = zeros(length(y),1);
n = 1;
w = zeros(Jb,1);

while n+block_size-1 <= length(y)
	Y = fft(y(n:n+block_size-1) .* window);
	w_in = GPSR_BB(Y, Bfft, 3, 'Verbose', 0, 'ToleranceA', 1);
	w = w_in*weight_b0 + w*weight_a1;
	y_re(n:n+block_size-1) = y_re(n:n+block_size-1) + (B'*w);
	
	n = n + hop_size;
end
