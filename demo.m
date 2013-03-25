% demo.m
% demo of modularized functions


addpath('bases');

[y fs] = wavread('tracks/al.wav');
y = y(:,1); % take first channel
y = y(1:fs*30); % take first few seconds

block_size = 1024;
hop_size = block_size/2;
window = hanning(block_size);

% oscillator basis
% B = sawtoothBasis(block_size, 200, fs, 1);

% wav file basis
% B = wavBasis('bases/spanish-phonemes');
% B = wavBasis('bases/alto-sax-free-jazz');
B = wavBasis('bases/piano2');
block_size = floor(size(B,2)/2);
hop_size = floor(block_size/16);
window = hanning(block_size);
% truncate basis vectors to block_size
B = B(:,1:block_size);

% pad out to hop size

tic
y_re = undercomplete2(y, B, block_size, hop_size, window, 0.0, 10);
toc

wavwrite(y_re, fs, 'constrained.wav');
