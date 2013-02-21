function [B] = sawtoothBasis(Nb, Jb, fs, detune)

if nargin < 4
	detune = 0;
end

% oscillators of varying frequency
B = zeros(Jb, Nb);
% construct basis by dividing 0-nyquist into equal parts
t = (0:Nb-1)/fs;
for jj=1:Jb
	freq = (jj-1)/Jb*fs/2 + detune;
	B(jj,:) = [1-2*rem(t*freq,1)];
end
