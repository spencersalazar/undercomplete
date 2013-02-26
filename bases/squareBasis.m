function [B] = squareBasis(Nb, Jb, fs, detune)

% 
% function [B] = squareBasis(Nb, Jb, fs, detune)
% 
% Generate Jb square-wave basis functions of length Nb with frequency from 0 to
% fs/2. 
% 
% === REQUIRED ARGUMENTS ===
% Nb
%	maximum signal length of any given basis function
% 
% Jb
%	number of basis functions to generate
% 
% fs
%   sampling rate
%
% === OPTIONAL ARGUMENTS ===
% detune
%	Value added to frequency of each basis function, to avoid comb-filter-like
%   coloration. 
%	(default: 0)
%
% === RETURN VALUE ===
% B
%	Jb*Nb matrix. Each row of B is a basis vector. 
%


if nargin < 4
	detune = 0;
end

% oscillators of varying frequency
B = zeros(Jb, Nb);
% construct basis by dividing 0-nyquist into equal parts
t = (0:Nb-1)/fs;
for jj=1:Jb
	freq = (jj-1)/Jb*fs/2 + detune;
	B(jj,:) = [1-2*round(rem(t*freq,1))];
end
