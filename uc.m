% uc.m

addpath('GPSR_6.0');

[y fs] = wavread('tracks/al.wav');
y = y(:,1); % take first channel


block_size = 1024;
hop_size = block_size/2;
window = hanning(block_size);

tic
% construct basis
% sawtooths of varying frequency
Nb = block_size; % basis length
Jb = 200; % number of basis functions
B = zeros(Jb, Nb);
% construct basis by dividing 0-nyquist into equal parts
for jj=1:Jb
	freq = (jj-1)/Jb*fs/2;
	t = (0:Nb-1)/fs;
	B(jj,:) = [1-2*fmod(t*freq,1)];
end
toc

% pad out to hop size
y = [y; zeros(hop_size - mod(length(y),hop_size),1)];
y_re = zeros(length(y),1);
n = 1;
tic
while n+block_size-1 <= length(y)
	w = GPSR_BB(y(n:n+block_size-1) .* window, B', 0, 'Verbose', 0, 'ToleranceA', 1);
	y_re(n:n+block_size-1) = y_re(n:n+block_size-1) + (B'*w);
	
	n = n + hop_size;
end
toc

wavwrite(y_re, fs, 'constrained.wav');
