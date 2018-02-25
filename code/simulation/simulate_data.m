function [x, y] = simulate_data(nobs, nrows, ncols)

mu1 = 1;
mu2 = 2;
sigma1 = 10;
sigma2 = 10;
class0size = floor(nobs/2);

x1 = randn([nrows, ncols, class0size])*sigma1 + mu1; %normrnd('norm', mu1, sigma1, [nrows, ncols, class0size]);
x2 = randn([nrows, ncols, nobs-class0size])*sigma2 + mu2; %normrnd('norm', mu2, sigma2, [nrows, ncols, nobs-class0size]);
x = cat(3, x1, x2);

y1 = zeros(class0size, 1);
y2 = ones(nobs-class0size, 1);
y = cat(1, y1, y2)+1; 

shuffled_order = randperm(nobs);
y = y(shuffled_order);
x = x(:,:,shuffled_order);


end