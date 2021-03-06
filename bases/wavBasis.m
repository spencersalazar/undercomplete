function [B] = wavBasis(directory, Nb)

%directory = 'bases/spanish-phonemes';

%
% function [B] = wavBasis(directory, Nb)
% 
%Nb = 0;
if nargin < 2
	Nb = -1;
end

paths = {};
max_size = 0;
sum_size = 0;
dirlist = dir(directory);
for ii = 1:length(dirlist)
	[pathstr filename ext] = fileparts(dirlist(ii).name);
	if strcmpi(ext, '.wav')
		thepath = fullfile(directory, dirlist(ii).name);
		%fprintf(1, 'path: %s\n', thepath);
		[samps chans] = wavread(thepath, 'size');
		if samps > max_size
			max_size = samps;
		end
		sum_size = sum_size + samps;
		paths{end+1} = thepath;
	end
end

if Nb == 0 % max size
	Nb = max_size;
end

if Nb == -1 % average size
	Nb = floor(sum_size / length(paths));
end

B = zeros(length(paths), Nb);

for ii = 1:length(paths)
%	fprintf(1, 'wavread path: %s\n', paths{ii});
	x = wavread(paths{ii});
	len = min(Nb,length(x));
	B(ii,1:len) = x(1:len)';
end
