%%--------------------------------------------------
%        Computing Transitions in the NGM
%%--------------------------------------------------

% This code computes the transition between two steady states after a
% once-and-for-all change in TFP. The framework is one of a Neoclassical
% Growth Model without endogenous labor supply. 

clear all

% Set Parameters

theta=0.67; % Labor Share
h=0.31; % Exogenous Labor
k_i=4; % Capital at initial steady state
i_i=0.25; % Investment at initial steady state
c_i=0.75; % Consumption at initial steady state
y_i = 1; % Normalized Output at intial steady state
delta = i_i/k_i; % Depreciation
z_ini = 1/((k_i.^(1 - theta) * h^(theta))^(1/theta)); % Labor productivity
beta = 1 / ((1-theta) * k_i^(-theta) * (z_ini*h)^(theta) + (1 - delta)); % Discount factor
T = 500; % Number of periods for transition

% Values for the final steady state.
k_f=8;
i_f=0.5;
c_f=1.5;
y_f=2;
z_fin=2*z_ini; %Shock


% Compute the transition

fun = @(x)euler(x, k_i, k_f, theta, delta, z_fin, h, beta, T);

x0 = ones(1,T).*(k_i);
options = optimoptions(@fsolve,'Algorithm','levenberg-marquardt',...
    'MaxFunctionEvaluations', 5000);
[k_values] = fsolve(fun, x0);
k_values(1)=k_i;
k_values(T)=k_f;

figure
plot(1:T, k_values(1:T))


function transition = euler(x, k_i, k_f, theta, delta, z_fin, h, beta, T)
    %k_0 = k_guess;
    %k_end = capital_ss;
    x(1) = k_i;
    x(T) = k_f;
    %transition = zeros(1,T);
    for j = 1:T-2
    %First period
        if j==1
            transition(j)=x(j+1)^(1-theta)*(z_fin*h)^(theta)+(1-delta)*x(j+1)-x(j+2)-...
                (k_i^(1-theta)*(z_fin*h)^(theta)+(1-delta)*k_i-x(j+1))*...
                beta*((1-theta)*x(j+1)^(-theta) * (z_fin*h)^(theta)+(1-delta));
        elseif j==(T-2)
            transition(j)=x(j+1)^(1-theta)*(z_fin*h)^(theta)+(1-delta)*x(j+1)-k_f-...
                (x(j)^(1-theta)*(z_fin*h)^(theta)+(1-delta)*x(j)-x(j+1))*...
                beta*((1-theta)*x(j+1)^(-theta) * (z_fin*h)^(theta)+(1-delta));
    %rest of periods:
        else
            transition(j)=x(j+1)^(1-theta)*(z_fin*h)^(theta)+(1-delta)*x(j+1)-x(j+2)-...
                (x(j)^(1-theta)*(z_fin*h)^(theta)+(1-delta)*x(j)-x(j+1))*...
                beta*((1-theta)*x(j+1)^(-theta) * (z_fin*h)^(theta)+(1-delta));
        end
    end
end