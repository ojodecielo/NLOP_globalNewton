% Author: Julia Schaeffer
% File name:  globnewton.m
% Globalized Newton method with Armijo stepsize rule
% It terminates when the norm of the gradient is below epsilon
% Function value <--- "func.m".
% Gradient value <--- "grad.m".
close all
% Input data
  x=[-1.2;1];     %initial point
  beta = .5;      %parameter beta for Armijo rule
  gamma = 1e-4;   %parameter gamma for Armijo rule
  alpha1 = 1e-6;  %alpha1 from globalized Newton method
  alpha2 = alpha1;%alpha2 from globalized Newton method
  p = .1;         %p from globalized Newton method
  epsilon=1e-9;   %stopping criterion
  obj=func(x);    %f(x)
  g=grad(x);      %grad(x)
  H=hesse(x);     %hesse(x)
  k=0;            %iteration step    
  s=g; 
  
fprintf('Initializing iterations: k = %g\n',k)
fprintf('Initial x = %s\n',sprintf('%f ', x))
fprintf('Initial f(x) = %f\n',obj)

% Begin method
  while  norm(g) > epsilon    %stopping criterion
    if det(H) == 0            %stop if Hessian singular
      break
    end  
    
    d = H\(-g);               %calculate d from H*d=-g
    norm_d = norm(d);         %euclidean norm of d
    norm_d_p = norm_d^p;
    norm_d_sq = norm_d^2;
    
    %Choice of s
    if -g'*d >= min(alpha1,alpha2*norm_d_p)*norm_d_sq
      s = d;
    else
      s = -g;
    end
    
    %Armijo step size rule
    sigma = 1;                %initial step size
    newobj = func(x + sigma.*s);
    while newobj > obj+gamma*sigma*g'*s     %searching for best sigma
      %fprintf('Trying sigma =%f\n',sigma)
      sigma = sigma*beta;
      newobj = func(x + sigma.*s);
    end
    
    %Updates
    x = x + sigma.*s;           %update x
    k = k + 1;                  %update iteration step
    fprintf('Iteration step k = %g\n',k)
    fprintf('Current x = %s\n',sprintf('%f ', x))
    obj=newobj;                 %update f(x)
    fprintf('Current f = %f\n',obj)
    g=grad(x);                  %update grad(x)
    H=hesse(x);                 %update Hessian(x)
  end


fprintf('Final iteration step k = %g\n',k)
fprintf('Minimizer x =  %s\n', sprintf('%f ', x))
fprintf('Minimal function value f(x) = %f\n',obj)
