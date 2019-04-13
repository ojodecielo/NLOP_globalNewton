%Author: Julia Schaeffer
% File name:  steepdesc.m
% Steepest descent method with Armijo stepsize rule
% It terminates when the norm of the gradient is below epsilon
% Function value <--- "func.m".
% Gradient value <--- "grad.m".

% Input data
  x=[-1.2;1];   %initial point
  beta = .5;    %parameter beta for Armijo rule
  gamma = 1e-4; %parameter gamma for Armijo rule
  obj=func(x);  %f(x)
  epsilon=1e-3; %stopping criterion
  g=grad(x);    %grad(x)
  k=0;          %iteration step                        

% Begin method
  while  norm(g) > epsilon    %stiooing criterion
    d = -g;                   %steepest descent direction
    sigma = 1;                %initial step size
    newobj = func(x + sigma.*d);
    
    %Armijo step size rule
    while newobj > obj+gamma*sigma*g'*d 
      sigma = sigma*beta;
      newobj = func(x + sigma.*d);
    end
    
    %if (mod(k,100)==1) fprintf('%5.0f %5.0f %12.5e \n',k,nf,obj); end
    
    %Updates
    x = x + sigma.*d;       %update x
    obj=newobj;             %update f(x)
    g=grad(x);              %update grad(x)
    k = k + 1;              %update iteration step
  end


fprintf('Iteration step k = %g\n',k)
fprintf('Minimizer x =  %s\n', sprintf('%f ', x))
fprintf('Minimal function value f(x) = %f\n',obj)
