clc
clear
format long;

rng(69)

%%%%%%%%%%%%%%%%%%%%
%%%% QUESTION 1 %%%%
%%%%%%%%%%%%%%%%%%%%

%Set Inputs for Black-Scholes
S0 = 499.75;    
K = 500;        
T = 1;          
r = 0.05;       
q = 0.01;
epsilon = 1;
sigma = 0.08;  
%sigma = 0.10; %Change for Question 8
%sigma = 0.06; %Change for Question 9

%Calculate Call Value with Black-Scholes Functions

d1 = getd1(S0, K, r, q, sigma, T);

d2 = getd2(d1, sigma, T);

V0 = BS(S0, K, r, q, T, d1, d2, epsilon);


fprintf('Premium (V0): %.4f\n', V0);

%%
%%%%%%%%%%%%%%%%%%%%
%%%% QUESTION 2 %%%%
%%%%%%%%%%%%%%%%%%%%

%Set Rebalancing Frequency to calculate Time Step dt using Function. For
%Question 6 simply change the Rebalancing Frecuency to the desired measurement.
%Also for Question 7 simply change mu to the desired value.

rebalancing_frequency = 'daily';
%rebalancing_frequency = 'weekly';
%rebalancing_frequency = 'monthly';
dt = getRebalancingFrequency(rebalancing_frequency);
M = 1000;
mu = 0.15;
%mu = -0.15; %Change for Question 7
%mu = 0; %Change for Question 7
Q = 50;

%Initialize and fill S matrix with all 1000 paths for the observed dollar
%and save ln of last values for each path.

S_paths = zeros(T/dt, M);
S_paths(1, :) = S0;
ln_ST = zeros(1, M);

for m = 1:M
    for i = 2:T/dt
        S_paths(i, m) = S_paths(i - 1, m) * exp((mu - sigma^2 / 2) * dt + sigma * sqrt(dt) * randn());
    end
    ln_ST(m) = log(S_paths(end, m));
end

Plot1 = figure;
histogram(ln_ST, Q, 'Normalization','pdf');

hold on;
x = linspace(min(ln_ST), max(ln_ST), T/dt);
y = normpdf(x, mean(ln_ST), std(ln_ST));
plot(x, y, 'r', 'LineWidth', 2);
hold off;

xlabel('Final Position ln(ST)');
ylabel('Probability Density');
title('Empirical Distribution vs. Theoretical Distribution');
legend('Empirical Distribution', 'Theoretical Distribution');
grid on;

%%
%%%%%%%%%%%%%%%%%%%%
%%%% QUESTION 3 %%%%
%%%%%%%%%%%%%%%%%%%%

%Initialize all relevant matrices
    
H_paths = zeros(T/dt, M);
B_paths = zeros(T/dt, M);
delta_paths = zeros(T/dt, M);
V_paths = zeros(T/dt, M);
Y_paths = zeros(T/dt, M);
    
%Set V_0 = H_0
    
V_paths(1, :) = V0;
H_paths(1,:) = V0;
    
%Calculate deltas with functions to get the time remaining to maturity and 
%Black-Scholes delta calculating function, and calculate B and H.
    
for m = 1:M
    for i = 1:T/dt
        T_remaining = getT_remaining(T, i, dt);
        delta_paths(i, m) = getBSdelta(S_paths(i, m), K, r, q, sigma, T_remaining, epsilon);
        B_paths(i, m) = H_paths(i, m) - delta_paths(i, m) * S_paths(i, m);
         if i < T/dt
             H_paths(i + 1, m) = delta_paths(i, m) * exp(q * dt) * S_paths(i + 1, m) + B_paths(i, m) * exp(r * dt);
         end
    end
end
    
%Calculate V with same functions as Question 1
    
for m = 1:M
    for i = 2:T/dt
        T_remaining = getT_remaining(T, i, dt);
        d1 = getd1(S_paths(i, m), K, r, q, sigma, T_remaining);
        d2 = getd2(d1, sigma, T_remaining);
        V_paths(i, m) = BS(S_paths(i, m), K, r, q, T_remaining, d1, d2, epsilon);
    end
end
    
%Calculate P&L Y
    
for m = 1:M
    for i = 1:T/dt
        if i < T/dt
            Y_paths(i, m) = (H_paths(i + 1, m) - H_paths(i, m)) - (V_paths(i + 1, m) - V_paths(i, m));
        end
    end
end


%%
%%%%%%%%%%%%%%%%%%%%
%%%% QUESTION 4 %%%%
%%%%%%%%%%%%%%%%%%%%

%Calculate net value X

X_paths = zeros(T/dt, M);

for m = 1:M
    for i = 1:T/dt
        X_paths(i, m) = H_paths(i, m) - V_paths(i, m);
    end
end


Plot2 = figure;
subplot(2, 1, 1);
plot(mean(Y_paths, 2), 'b', 'LineWidth', 2);
hold on;
plot(std(Y_paths, 0, 2), 'r', 'LineWidth', 2);
hold off;
xlabel('Time (Trading Days)');
%xlabel('Time (Weeks)'); %Change for Question 6
%xlabel('Time (Months)'); %Change for Question 6
ylabel('Value');
title('Mean and Standard Deviation P&L Y over Time');
legend('Mean', 'Standard Deviation');
grid on;

subplot(2, 1, 2);
plot(mean(X_paths, 2), 'b', 'LineWidth', 2);
hold on;
plot(std(X_paths, 0, 2), 'r', 'LineWidth', 2);
hold off;
xlabel('Time (Trading Days)');
%xlabel('Time (Weeks)'); %Change for Question 6
%xlabel('Time (Months)'); %Change for Question 6
ylabel('Value');
title('Mean and Standard Deviation Net Value X over Time');
legend('Mean', 'Standard Deviation');
grid on;

%%
%%%%%%%%%%%%%%%%%%%%
%%%% QUESTION 5 %%%%
%%%%%%%%%%%%%%%%%%%%

%Get deltas for each relevant rebalancing date

RebalancingDay_125 = delta_paths(125,:);
RebalancingDay_252 = delta_paths(end,:);

Plot3 = figure;

subplot(2, 1, 1);
histogram(RebalancingDay_125, Q, 'Normalization', 'pdf');
xlabel('Delta at 125th Rebalancing Day');
ylabel('Probability Density');
title('Distribution of Delta at the 125th Rebalancing Day');
grid on;

subplot(2, 1, 2);
histogram(RebalancingDay_252, Q, 'Normalization', 'pdf');
xlabel('Delta at Maturity');
ylabel('Probability Density');
title('Distribution of Delta at Maturity');
grid on;



%%
%%%%%%%%%%%%%%%%%%%%%
%%%% QUESTION 10 %%%%
%%%%%%%%%%%%%%%%%%%%%

%Set Inputs for Black-Scholes
S0 = 499.75;    
K = 500;        
T = 1;          
r = 0.05;       
q = 0.01;       
sigma = 0.10; %Set market volatility  
epsilon = 1;

%Calculate Call Value with Black-Scholes Functions

d1 = getd1(S0, K, r, q, sigma, T);

d2 = getd2(d1, sigma, T);

V0 = BS(S0, K, r, q, T, d1, d2, epsilon);


fprintf('Premium (V0): %.4f\n', V0);

rebalancing_frequency = 'daily';
dt = getRebalancingFrequency(rebalancing_frequency);
M = 1000;
mu = 0.15;
Q = 50;

S_paths = zeros(T/dt, M);
S_paths(1, :) = S0;
ln_ST = zeros(1, M);

for m = 1:M
    for i = 2:T/dt
        S_paths(i, m) = S_paths(i - 1, m) * exp((mu - sigma^2 / 2) * dt + sigma * sqrt(dt) * randn());
    end
    ln_ST(m) = log(S_paths(end, m));
end
    
H_paths = zeros(T/dt, M);
B_paths = zeros(T/dt, M);
delta_paths = zeros(T/dt, M);
V_paths = zeros(T/dt, M);
Y_paths = zeros(T/dt, M);
    
V_paths(1, :) = V0;
H_paths(1,:) = V0;
    
for m = 1:M
    for i = 2:T/dt
        T_remaining = getT_remaining(T, i, dt);
        d1 = getd1(S_paths(i, m), K, r, q, sigma, T_remaining);
        d2 = getd2(d1, sigma, T_remaining);
        V_paths(i, m) = BS(S_paths(i, m), K, r, q, T_remaining, d1, d2, epsilon);
    end
end
    
sigma = 0.08; %Set trader's estimated volatility
    
for m = 1:M
    for i = 1:T/dt
        T_remaining = getT_remaining(T, i, dt);
        delta_paths(i, m) = getBSdelta(S_paths(i, m), K, r, q, sigma, T_remaining, epsilon);
        B_paths(i, m) = H_paths(i, m) - delta_paths(i, m) * S_paths(i, m);
         if i < T/dt
             H_paths(i + 1, m) = delta_paths(i, m) * exp(q * dt) * S_paths(i + 1, m) + B_paths(i, m) * exp(r * dt);
         end
    end
end

%Calculate P&L Y
    
for m = 1:M
    for i = 1:T/dt
        if i < T/dt
            Y_paths(i, m) = (H_paths(i + 1, m) - H_paths(i, m)) - (V_paths(i + 1, m) - V_paths(i, m));
        end
    end
end

%Calculate net value X

X_paths = zeros(T/dt, M);

for m = 1:M
    for i = 1:T/dt
        X_paths(i, m) = H_paths(i, m) - V_paths(i, m);
    end
end

Plot4 = figure;
subplot(2, 1, 1);
plot(mean(Y_paths, 2), 'b', 'LineWidth', 2);
hold on;
plot(std(Y_paths, 0, 2), 'r', 'LineWidth', 2);
hold off;
xlabel('Time (Trading Days)');
ylabel('Value');
title('Mean and Standard Deviation P&L Y over Time');
legend('Mean', 'Standard Deviation');
grid on;

subplot(2, 1, 2);
plot(mean(X_paths, 2), 'b', 'LineWidth', 2);
hold on;
plot(std(X_paths, 0, 2), 'r', 'LineWidth', 2);
hold off;
xlabel('Time (Trading Days)');
ylabel('Value');
title('Mean and Standard Deviation Net Value X over Time');
legend('Mean', 'Standard Deviation');
grid on;

