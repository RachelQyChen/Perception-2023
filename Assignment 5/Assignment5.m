%% Perception 2023 Assignment 5

%% 2 a)
clear; close all; clc;
tRadius= 1; % radius of the target
x = randi(10); % The x-coordinate of the center of the target
y = randi(10); % The y-coordinate of the center of the target
trials = 1e6; 

var1 = 1; % variance of both Gaussians for cue1
sig1 = sqrt(var1);
var2 = 4; % variance of both Gaussians for cue2
sig2 = sqrt(var2);



% Drawing samples from Gaussians
X1 = x+ mvnrnd(0, sig1, trials);
X2 = x+ mvnrnd(0, sig2, trials);
Y1 = y+mvnrnd(0, sig1, trials);
Y2 = y +mvnrnd(0, sig2, trials);

%% a)
% If Michael uses only cue 1, then the samples drawn from the Gaussian X1
% and Y1 are the x and y coordinates of the estimates of the center of the
% target. The estimate is correct if it lies inside the circle of radius 1
% around the center. The radial distance of each estimate from the center
% can be computed using the Euclidean distance formula:
%%
% $$d = \sqrt{(x - \bar{x})^2 + (y - \bar{y})^2} $$
%%
% And the estimate is correct if this radial distance is less than or equal
% to the radius of the circle which is 1.
%%
% The probability of correct estimate can then be computed as the ratio of
% the count of the correct estimates by the total number of trials.
%%
figure();
scatter(X1, Y1, 1)
hold on;
plot(x, y, 'r*')
xlabel('x coordinations')
ylabel('y coordinations')
title('Scatter plot of estimates with only cue_1')
axis equal

distance = sqrt((X1 - x).^2 + (Y1 - y).^2);
cCorrect = sum(distance <= tRadius);
pCorrect = cCorrect/trials;

figure();
scatter(X2, Y2, 1)
hold on;
plot(x, y, 'g*')
xlabel('x coordinations')
ylabel('y coordinations')
title('Scatter plot of estimates with only cue_2')
axis equal
distance2 = sqrt((X2 - x).^2 + (Y2 - y).^2);
cCorrect2 = sum(distance2 <= tRadius);
pCorrect2 = cCorrect2/trials;

%% 
w = var2/(var1 + var2);
v = var2/(var1 + var2);

Xb = w * X1 + (1 - w) * X2;
Yb = v .* Y1 + (1 - v) .* Y2;

figure();
scatter(Xb,Yb, 1)
hold on;
plot(x, y, 'r*')
xlabel('x coordinations')
ylabel('y coordinations')
title('Scatter plot of estimates with combined cue')
axis equal

distance44 = sqrt((Xb - x).^2 + (Yb - y).^2);
cCorrect44 = sum(distance44 <= tRadius);
pCorrect44 = cCorrect44/trials;


%% Question 3a
x_all= 1:100;
Denis_prior_mean = 50; 
Denis_prior_sd = sqrt(64); 

prior_unnormalized = normpdf(x_all, Denis_prior_mean, Denis_prior_sd);
%prior_normalized_1 = prior_unnormalized ./ max(prior_unnormalized);
prior_normalized_2 = prior_unnormalized ./ sum(prior_unnormalized);

likelihood_mean = 30;
likelihood_sd = sqrt(100);
likelihood_unnormalized = normpdf(x_all, likelihood_mean, likelihood_sd);
%likelihood_normalized_1 = likelihood_unnormalized ./ max(likelihood_unnormalized);
likelihood_normalized_2 = likelihood_unnormalized ./ sum(likelihood_unnormalized);

posterior_unnormalized = prior_normalized_2 .* likelihood_normalized_2;
%posterior_unnormalized = prior_normalized_1 .* likelihood_normalized_1;
posterior_normalized_2 = posterior_unnormalized ./ sum(posterior_unnormalized);

figure();
plot(x_all, prior_normalized_2, '-.','DisplayName', 'prior', 'LineWidth', 2)
hold on;
plot(x_all, likelihood_normalized_2, ':', 'DisplayName', 'likelihood', 'LineWidth', 2)
plot(x_all, posterior_normalized_2, 'DisplayName', 'posterior', 'LineWidth', 2)
xlabel('All possible points')
ylabel('Probability')
ylim([0, 0.1])
legend()

%% Question 3b
x_all= 1:100;
Denis_prior_mean = 50; 
Denis_prior_sd = sqrt(64); 

prior_unnormalized = normpdf(x_all, Denis_prior_mean, Denis_prior_sd);
%prior_normalized_1 = prior_unnormalized ./ max(prior_unnormalized);
prior_normalized_2 = prior_unnormalized ./ sum(prior_unnormalized);

likelihood_mean = 30;
likelihood_sd = sqrt(36);
likelihood_unnormalized = normpdf(x_all, likelihood_mean, likelihood_sd);
%likelihood_normalized_1 = likelihood_unnormalized ./ max(likelihood_unnormalized);
likelihood_normalized_2 = likelihood_unnormalized ./ sum(likelihood_unnormalized);

posterior_unnormalized = prior_normalized_2 .* likelihood_normalized_2;
posterior_normalized_2 = posterior_unnormalized ./ sum(posterior_unnormalized);

figure();
plot(x_all, prior_normalized_2, '-.','DisplayName', 'prior', 'LineWidth', 2)
hold on;
plot(x_all, likelihood_normalized_2, ':', 'DisplayName', 'likelihood', 'LineWidth', 2)
plot(x_all, posterior_normalized_2, 'DisplayName', 'posterior', 'LineWidth', 2)
xlabel('All possible points')
ylabel('Probability')
ylim([0, 0.1])
legend()

%% Question 4

clear; close all; clc;

trials = 1e6; 
x = randi(10).*ones(trials, 1); % x center of all targets
y = randi(10).*ones(trials, 1); % y center of all targets
tRadius= ones(trials, 1); % All radius for ragets
possible_cues = 1:100; 
var = 10.*ones(trials, 1); % All variance of targets
%sig = sqrt(10).*ones(trials, 1);
sig = sqrt(10);

money = 1000;
cost = 10;


%empty container
x_choice = nan(trials, length(possible_cues));
y_choice = nan(trials, length(possible_cues));


X_cue = x + sig.*randn(trials, length(possible_cues));
Y_cue = y + sig.*randn(trials, length(possible_cues));


%Generate weighted results
for jj = possible_cues
    w = ones(jj, 1)./jj;
    v = ones(jj, 1)./jj;
    x_choice(:,jj) = X_cue(:, 1:jj)*w;
    y_choice(:,jj) = Y_cue(:, 1:jj)*w;
end 

distance_3b = sqrt((x_choice - x).^2 + (y_choice - y).^2);
cCorrect3 = sum(distance_3b <= tRadius);
pCorrect3 = cCorrect3/trials;

gain = 1000.*pCorrect3; 
lose = cost*possible_cues;

money_in_pocket = gain - lose;
[max_money_gain, max_count] = max(money_in_pocket);

figure();
plot (possible_cues,money_in_pocket, 'LineWidth', 1.2, 'Color', 'r'); hold on;
plot (possible_cues(max_count), max_money_gain, 'o', 'MarkerSize', 8, 'MarkerFaceColor','r')
text(possible_cues(max_count)-10, max_money_gain+30, 'critical value of cue: 32', 'Color', 'r')
xlabel('Number of possible cues')
ylabel('Money earning')
ylim([-100 600]);
yline(0, ':', 'Color', 'b');
