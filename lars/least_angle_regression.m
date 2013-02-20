% least angle regression
% @author: Sam Abeyruwan

function least_angle_regression(Data)

% Standardize the data to have zero mean and unit L2-norm
% and the response is centered around 0; so that the bias term in not
% needed

close all;
[m n] = size(Data)

X = zeros(m, n - 1);
y = zeros(m, 1);
%
for j = 1 : n - 1
	% remove the mean
	X(:, j) = Data(:, j) - mean(Data(:, j));
	% unit L2 - norm
	div = sqrt(X(:, j)' * X(:, j));
	X(:, j) = X(:, j) / div;
end
%
y = Data(:, n) - mean(Data(:, n));
%

printf("Forward Stagewise Linear Regression \n");


epsilon = 0.8;
cnt = 1;
stagewise_steps = 2000;

% candiate vector of regression coefficientes
beta_hat = zeros(n - 1, 1);
% regression indexes
j_hat = [];

r = y;
while (stagewise_steps >= cnt)
	curr_cor = X' * r;
	[_, idx] = max(abs(curr_cor));
	beta_hat(idx) = beta_hat(idx) + epsilon * sign(curr_cor(idx));
	r = r - epsilon * sign(curr_cor(idx)) * X(:, idx);	
	if (sum(ismember(j_hat, idx)) == 0), j_hat = [j_hat, idx]; end	
	cnt = cnt + 1;	
end
%disp(j_hat)
%disp(beta_hat(find(beta_hat)));

y_hat = X * beta_hat;
errorv = sqrt((y - y_hat)' * (y - y_hat));
printf("error1=%f \n", errorv);
f1 = figure;
hold on;
plot([1:100], y([1:100]), '-y',[1:100], y_hat([1:100]), '-or');

beta2 = inv(0.01* eye(n - 1, n -1) + X' * X) * X' * y;
y_hat2 = X * beta2;
errorv = sqrt((y - y_hat2)' * (y - y_hat2));
printf("error2=%f \n", errorv);

plot([1:100], y_hat2([1:100]), '-xg'); 


printf("lars OLS\n");
[beta3, A] = lars(Data(:, 1:n - 1), Data(:, n), 'lars', 1500)
y_hat3 = X * beta3;
errorv = sqrt((y - y_hat3)' * (y - y_hat3));
printf("error3=%f \n", errorv);
plot([1:100], y_hat3([1:100]), '-xb'); 
legend show Location NorthEastOutside;
title("A sample run");
legend("True response", "ISR t=Inf", "L2-norm lambda=0.01", "LARS t=1500");

print -djpg /home/sam/tmp/octave_test2/lars/plot1.jpg

hold off;




beta_lars = [];
beta_lasso = [];

inc = 1;
max_v = 3600;
len = length([0 :100:max_v]);
for t = 0 : 100: max_v
	[curr_beta, A] = lars(Data(:, 1:n - 1), Data(:, n), 'lars', t)
	[curr_beta22, A22] = lars(Data(:, 1:n - 1), Data(:, n), 'lasso', t)
	beta_lars(inc, :) = curr_beta';
	beta_lasso(inc, :) = curr_beta22';
	inc = inc + 1;	
end

[mm nn] = size(beta_lars)
% disp(beta_lars);

NumberOfPlots= nn;
ColorSet=varycolor(NumberOfPlots);

f2  = figure;
hold on;
for r = 1 : nn	
	plot([0 :100:max_v]', beta_lars(:, r), 'linestyle', "-.", 'Color',ColorSet(r, :));	
end
legend('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', 'Location', 'SouthWestOutside');
title('LARS traces');
print -djpg /home/sam/tmp/octave_test2/lars/plot2.jpg
hold off;

f3 = figure;
hold on;
for r = 1 : nn	
	plot([0 :100:max_v]', beta_lasso(:, r), 'Color',ColorSet(r, :));	
end
legend('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', 'Location', 'SouthWestOutside');
title('LASSO traces');
print -djpg /home/sam/tmp/octave_test2/lars/plot3.jpg
hold off;

