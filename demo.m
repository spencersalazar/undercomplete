% demo.m
% demo of modularized functions


addpath('bases');

[y fs] = wavread('tracks/bowie.wav');
y = y(:,1); % take first channel

block_size = 1024;
hop_size = block_size/2;
window = hanning(block_size);

B = sawtoothBasis(block_size, 200, fs, 1);

% pad out to hop size

y_re = undercomplete(y, B, block_size, hop_size, window);

wavwrite(y_re, fs, 'constrained.wav');
