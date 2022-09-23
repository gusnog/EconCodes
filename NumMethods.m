%% Simple Numerical Methods. Based on Fehr and Kindermann 
% Root finding methods. 

% Normally we can solve a linear system by applying simple linear algebra
% methods. However, when the system is not linear and/or we do not have an
% analytical solution for the derivative of the function we are trying to
% find the solution for, we can apply some numerical techniques to find the
% solution to our problem. 

%This introduces some methods in the simplest one-dimensional case: 

%% Bisection Algorithm in One Dimension 

%This method takes to reference points in the domain of a function and
%iteratively tries to find a root (solution) for it. We will start with two
%point, a and b. Then, we define a point x which will be in the middle
%point of the segment between these two. Then, depending on the side of the
%function we are in, we will update our algorithm in one direction or the
%other. 

% Define the parameters and functions to use: 

a = 0.05; 
b = 0.25;

% Define the value of the function evaluated in each one of our points:
fa = 0.5*(a)^(-0.2)+0.5*(a)^(-0.5)-2;
fb = 0.5*(b)^(-0.2)+0.5*(b)^(-0.5)-2;

%Now we check if in the interval [a,b] we have one root. If not, we need to
%redefine our boundary points. How we check this? Think that we have a
%decreasing function, the root is at the point where the function becomes
%(numerically) zero. This is at the point where it crosses the x-axis.
%Then, if the function evaluated in a and in b changes signs
%this means that at some point between [a,b] we have the crossing point. If
%the product of the function evaluated at a and evaluated at b was of
%positive sign, this means that a and b are point which image is on the
%positive side of the function (above the x-axis), or both points have
%associated images below the x-axis. This cannot be the case if we want to
%find the root by bisecting the interval between these points. 

%Check this: 

if (fa*fb >= 0)
    disp('There is no root in [a,b]')
else
    disp('The root is in [a,b]')
end

%Then, once checked we can go to our algorithm:

%Define the counter for iterations, together with the tolerance levels 
% and the maximum number of iterations 
iter=0; 
tol = 1e-6; 
maxit = 200; 

%Our iterative algorithm will keep updating the reference value for
%bisection until we eithe rreach convergence or exhaust the number of
%iterations. After each iteration, the algorithm will check if we are
%sufficiently close to a solution. This will be done by comparing the
%current iteration with the past one. If the difference between them is
%sufficiently small (1e-6), then we stop the algorithm. Otherwise, we
%update our bisection value by taking new values for a and b. 

while (iter < maxit)
    iter=iter+1;
%First calculate the middle point of the segment [a,b]. This will be the
%point at which we will evaluate the function:
x = (a+b)/2;

%Notice that each iteration of the algorithm will be halving the interval.

%Then, evaluate our function at that value:
fx= 0.5*(x)^(-0.2)+0.5*(x)^(-0.5)-2;

%We check convergence:
if (abs(x-a) < tol)
    disp(['Convergence achieved at iteration: ' num2str(iter)])
    disp(['The value of the root is: ' num2str(x)])
    disp(['The value of the function at the root is: ' num2str(fx)]) 
break    
else  %if no convergence, then update our values:
    %If a and our bisection value are on different sides of the function
    %(taking the root as reference), then the root must be closer to b (the
    %other reference value), since we have checked that a and b are at
    %different sides of the function taking the root as reference. Then,
    %we make our reference point equal to b, and repeat the procedure.
    if (fa*fx < 0)
        b=x;
        fb=fx;
    else %If fa*fx is positive, this means that by taking a as a reference we will be moving
        %closer to the root coming from a's side. Then, update the value
        %taking a as a reference. It is all about knowing form what side
        %should we approach our bisection point. 
        a=x;
        fa=fx;
    end
end
end


%% Newton's Method in One Dimension:

% The idea of Newton's method is to approximate the solution to our problem
% through successive linearizations of the corresponding function. This
% allows to change a nonlinear root finding problem into a sequence of
% linear problems. The departing point is working with a first order Taylor
% approximation of the function we want to solve for. 

%The main idea of this is use the gradient of the function to approximate
%the solution. This requires taking a reference value x_0. Then get the
%value of the function at this point, and evaluate its derivative at this
%point as well. This allows us to get a new reference value, x_1, that
%should be closer to the root of the function. We repeat this approximation
%until we reach some tolerance value: 

%First define our reference value: 

xold = 0.05; 

%Then initialize the iteration counter: 
iter=0; 
maxit = 200; 
tol = 1e-6;

%Initialize our loop: 

while (iter < maxit)
 %First calculate the value of the function at the initial reference point:
iter=iter+1;
 f= 0.5*(xold)^(-0.2)+0.5*(xold)^(-0.5)-2;

 %We can compute analitically the derivative by hand. Having this we can
 %then write it and calculate the value of the derivative of the function
 %at our reference point:

 fprime = -0.1*(xold)^(-1.2)-0.25*(xold)^(-1.5)-2;

 %Then update our reference value:
 x = xold - (f/fprime);

 %Now, check convergence:
 if (abs(x-xold) < tol) %Convergence achieved if the updated value is 
     % sufficiently close to the reference value we started the iteration
     % with
     disp(['Convergence achieved after iteration: ' num2str(iter)])
     disp(['The value of the root is: ' num2str(x)])
     break
 else %If still no convergence, then take our updated value as the new 
     % reference for starting the iteration:
     xold=x;
 end
end

%% Secant Method

%This is inspired in Newton's Algorithm. Used in situations in which
%computing analitically the derivative of the function we are interested in
%is difficult or unfeasible. In this case, we approximate the derivative by
%a secant. 

%The secant will be such that:
% fprime = (f(x_i)-f(x_{i-1}))/(x_i-x_{i-1})

x_old = 0.05;
x_curr = 0.25;

iter=0; 
maxit = 200; 
tol = 1e-6;

while (iter < maxit)
iter=iter+1;
    f1 = 0.5*(x_old)^(-0.2)+0.5*(x_old)^(-0.5)-2;
    f2 = 0.5*(x_curr)^(-0.2)+0.5*(x_curr)^(-0.5)-2;
  
    %Then update our reference value:
 x = x_curr - ((x_curr-x_old)/(f2-f1))*f2 ;

  if (abs(x-x_curr) < tol) %Convergence achieved if the updated value is 
     % sufficiently close to the reference value we started the iteration
     % with
     disp(['Convergence achieved after iteration: ' num2str(iter)])
     disp(['The value of the root is: ' num2str(x)])
     break
  else %If still no convergence, then take our updated value using the secant 
      % and then use it as our current reference. 
     x_old=x_curr;
     x_curr=x;
 end
end

%Generally, the secant method is more efficient than Newton's algorithm.

%% Fixed Point Iteration Methods

%This method exploits the fact that we can write the for f(x)=0:
% x=x+sigma*f(x)==g(x). 
%That is, use the fact that we have, for some x, that f(x)=0. Then, we can
%multiply by a constant sigma at both sides, and add x to write a new
%funtion g(x)=x+sigma*f(x). This function, provided that sigma is not zero,
%will mean that the solution x* is a fixed point to this equation.

%Then, we are now updating according to this factor sigma. Generally, the
%lower this value of sigma is, the slower would the convergence be.
%However, it can be the case that sigma is too large and our method does
%not converge at all. Try sigma 0.5 and likely it will not.

%Naturally, other conditions must be satisfied in order to achieve
%convergence. Particularly, if the derivative of g(x) evaluated at the
%solution is smaller than one, and if the initial guess is sufficiently
%close to the root of f(x), convergence will be achieved. 

%Proceed similarly as before. Just adding the sigma factor to update our x:

x_old=0.12;
sigma = 0.2;

iter=0; 
maxit = 200; 
tol = 1e-6;

while (iter<maxit)
  
iter=iter+1;
    %First evaluate our function on the initial value:
f = 0.5*(x_old)^(-0.2)+0.5*(x_old)^(-0.5)-2;

%Update the value of x:

x = x_old+sigma*f;

%Check convergence:
  if (abs(x-x_curr) < tol) %Convergence achieved if the updated value is 
     % sufficiently close to the reference value we started the iteration
     % with
     disp(['Convergence achieved after iteration: ' num2str(iter)])
     disp(['The value of the root is: ' num2str(x)])
     break
  else %If still no convergence, then take our updated value  
      % and then use it as our current reference for the next iteration: 
     x_old=x;
  end


end

if iter >= maxit
   disp('Convergence not achieved');
end


