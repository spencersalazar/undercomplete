function [beta, A] = lars(X, y, option, t)

% Least Angle Regression (LARS) algorithm
% @author: Sam Abeyruwan

if nargin < 4, t = Inf; end
if nargin < 3, option = 'lars'; end
if strcmpi(option, 'lasso'), lasso = true; else, lasso = false; end
	
printf("running mode=%s  t=%f \n", option, t);	
	
eps = 1e-10; % effective zero

[m, n] = size(X)
% preprocessing
% covariates needs to be standardized to have mean 0 and unit L2-norm
% and is centered around 0 (response has mean 0)
% This is needed because the the project is obtained directly from the dot product
for j = 1 : n
	X(:, j) = X(:, j) - mean(X(:, j));
	X(:, j) = X(:, j) / sqrt(X(:, j)' * X(:, j));
end	
y = y - mean(y);

sup_active = min(n, m - 1); % maximum number of variables in the final active set
len_t = length(t);

beta = zeros(n, 1); % candidate regression coefficients
mu = zeros(m, 1); % predication vector

valid_sign = true;
active = 0; % variable that counts the max rounds
initialize = true;

A = []; % subset of indices of the active set
sign_A = []; % sign stays constant as the active set increases

c = [];
sup_c = 0;
sup_c_idx = 0;

% LARS steps length
gamma_active = 0;


% restrictions
stop_loop = false;
norm1_curr = 0;
norm1_prev = 0;

% LARS loop
while active < sup_active, 
	% Eq. 2.8
	c = X' * (y - mu); % current correlation
	sup_c =  max(abs(c)); % greatest absolute current correlation
	
	if sup_c < eps || isempty(t), break; end % early stopping criteria
	% intialize the activity
	if initialize, 
		[_, sup_c_idx] = max(abs(c)); 
		initialize = false;
	end	
	
	if valid_sign, 
		A = [A, sup_c_idx]; % add the covariate index to active set
		active = active + 1;
	end	
	
	% Eq. 2.10
	sign_A = sign(c(A)); % sign {-1, 1}
	
	complement_A = setdiff(1:sup_active, A); % complement set or inactive set
	complement_A_size = length(complement_A);
	
	% Eq. 2.4
	X_A = X(:, A);
	% Eq. 2.10
	for j = 1 : length(A)
		X_A(:, j) = sign_A(j) * X_A(:, j);	
	end
	% Eq. 2.5
	G_A = X_A' * X_A; % Gram matrix
	I_A = ones(length(A), 1);
	INV_G_A = inv(G_A);
	A_A = 1 / sqrt(I_A' * INV_G_A * I_A); 
	% Eq. 2.6
	w_A = A_A * INV_G_A * I_A;
	u_A = X_A * w_A; % equianglular vector
	% Eq. 2.11
	a = X' * u_A; % inner product vector	
		
	if active == sup_active,
		gamma_active = sup_c / A_A; % exceptional last stage
	else
		gamma_hat = zeros(complement_A_size, 2); % min of two
		for  j_itr = 1 : complement_A_size
			j = complement_A(j_itr);
			gamma_hat(j, :) = [( (sup_c - c(j)) / (A_A - a(j)) ), ( (sup_c + c(j)) / (A_A + a(j)) )];	
		end
		[gamma_active, sup_c_idx, _] = minplus(gamma_hat);
	end	
	
	beta_tmp(A) = beta(A) + gamma_active * diag(sign_A) * w_A; % update the coefficient estimate	
	stop_loop = false;	
	A_old = A;
	if strcmpi(option, 'lasso'),
	        % The LARS-Lasso relationship
		valid_sign = true; 
		d_hat = diag(sign_A' * w_A);
		% Eq. 3.4
		gamma_cross = -beta(A) ./ d_hat;		
		% Eq. 3.5
		[gamma_wave, gamma_wave_min_idx,  _] = minplus(gamma_cross);
		% Lasso modification (Theorem 1)
		if isnan(gamma_wave),
			gamma_wave = Inf;	
		end
		
		if gamma_wave < gamma_active,
			% Eq. 3.6			
			gamma_active = gamma_wave;
			beta_tmp(A) = beta(A) + gamma_active * diag(sign_A) * w_A; % Lasso update
			A(gamma_wave_min_idx) = []; % delete the zero-crossing element
			active = active - 1;
			valid_sign = false;								
		end	
	end		
	
	if t ~= Inf,
		norm1_curr = norm(beta_tmp(A_old), 1);
		if norm1_prev <= t && norm1_curr >= t,			
			% Eq. 5.17			
			beta_tmp(A_old) = beta(A_old) + A_A * (t - norm1_prev) * diag(sign_A) * w_A;			
			stop_loop = true;	                                                                               		
		end
		norm1_prev = norm1_curr;
	end
	
	mu = mu + gamma_active * u_A; % update the predication vector			
	beta(A_old) = beta_tmp(A_old);
	if stop_loop, break; end			
end	

	