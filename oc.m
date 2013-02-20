% oc.m

addpath('GPSR_6.0');
addpath('lasso');
addpath('lasso/sub');

[y fs] = wavread('ohno.wav');

% construct basis
% sawtooths of varying frequency
Nb = 1024; % basis length
Jb = 100; % number of basis functions
B = zeros(Jb, Nb);
% construct basis by dividing 0-nyquist into equal parts
for jj=1:Jb
	freq = (jj-1)/Jb*fs/2;
	t = (0:Nb-1)/fs;
	B(jj,:) = [1-2*fmod(t*freq,1)];
end

y = y(1:Nb,1); % take first Nb samples of first channel

%[w A] = lars(B,y,'lasso');
% 
% [w resnorm residual exitflag] = lsqnonneg(B',y');
% for ii = 1:1000
% 	[w resnorm residual exitflag] = lsqnonneg(B',y', w);
% end

w = GPSR_BB(y,B',1);

y_re = B'*w;

wavwrite(y_re, fs, 'constrained.wav');
wavwrite(y, fs, 'orig.wav');
