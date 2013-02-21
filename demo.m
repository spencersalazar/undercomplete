% demo.m
% demo of modularized functions


addpath('bases');

[y fs] = wavread('tracks/al.wav');
y = y(:,1); % take first channel

block_size = 1024;
hop_size = block_size/2;
window = hanning(block_size);

B = sawtoothBasis(block_size, 200, fs, 1);

% pad out to hop size

tic
y_re = undercomplete2(y, B, block_size, hop_size, window, 0.8);
toc

wavwrite(y_re, fs, 'constrained.wav');
