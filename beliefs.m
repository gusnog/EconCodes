%% The effect of observation noise in social learning (Based on Chamley's Rational Herds. Ch.3, 2004)
% Simple illustration of what happens when an agent's actions do not 
% convey all the information about his private information/beliefs. 
% In the benchmark, with the actions of the individual perfectly revealing
% his private signal, we have that the precision of the public belief
% (which is updated every period after an agent takes his actions)
% increases faster than when there is observation noise about the private
% signal an agent that moved in that period received.


% Basic timing: 
% 1) Given some initial public beliefs about the aggregate state, 
% nature draws a value for an aggregate parameter of interest (theta,
% for example)
% 2) An agent receives a private signal, s. Updates his private beliefs,
% using the public belief as his prior. 
% 3) That agent makes a decision to maximize his expected payoff
% 4) All agents, upon seeing agent t's actions, update the public belief.
% This is what we care about. How these beliefs will evolve, as in the next
% period an agent will update his private beliefs based on this and his own
% private signal. 

% Note: The aggregate state random variable, the public belief, and the
% private signal are all Gaussian. 

%% Without Information Noise --> An agent's actions perfectly reveal his private signal.

%Create a vector rho:
rho=zeros(100,1);

rho_theta = 1;
rho_epsi = 1;

for i=1:100
  rho(i) = rho_theta + (i-1)*rho_epsi;
end

%If we were to write it in terms of the variances instead of the
%precisions:

sigma=zeros(100,1);
sigma_theta = 1/rho_theta;
sigma_epsi = 1/rho_epsi;


for i=1:100
  sigma(i) = sigma_theta + (i-1)*sigma_epsi;
end



%% With observational noise --> Action of an agent does not fully reveal his private signal.

%The public beliefs now evolve according to:
% rho_t+1 = rho_t + 1/(sigma_epsilon^2 + sigma_eta^2*(1+rho_t+sigma_epsilon^2)^2)

% We define the values for the parameters: 

sigma_theta_noise = 1;
sigma_epsi_noise = 0.75; 
sigma_eta_noise = 0.423;

% Create the vector rho, for the public belief

rho_noise = zeros(100,1);

rho_noise(1) = 1; %Initial value, equal to rho_theta. 

% Updating Rule

for t=2:100
    rho_noise(t) = rho_noise(t-1) + 1/(sigma_epsi_noise^2 + sigma_eta_noise^2*(1+rho_noise(t-1)*sigma_epsi_noise^2)^2);
end    

%% Conclusion

% Then, without information noise, the public belief increases by a
% constant amount rho_epsi every period. It will be a linear function of
% the number of observations. With observation noise, as the number of
% periods goes to infinity the increment in the precision across periods
% tends to zero (it will be a concave function of the number of periods).
% To check this visually, see the following graph.

t = linspace(1,100);

figure;
pl=plot(t, rho, t, rho_noise), axis tight; title('Evolution of Public Beliefs'); xlabel('Observations'), ylabel('Precision');
set(pl,'LineWidth',2);
legend({'Without Observation Noise','With Observation Noise'},'Location','best')