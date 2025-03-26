N = 1000;
mu = [0 0 0];
sigma = [1 0 0;
        0 1 0;
        0 0 1];

count12 = 0;
count13 = 0;
for i = 1:N
    R = mvnrnd(mu, sigma, 100);
    H12 = ttest(R(:, 1), R(:, 2));
    if H12
        count12 = count12 + 1;
    end
    H13 = ttest(R(:, 1), R(:, 3));
    if H13
        count13 = count13 + 1;
    end
end

prob12 = count12/N;
prob13 = count13/N;
