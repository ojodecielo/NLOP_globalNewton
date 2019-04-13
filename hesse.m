%File name:  hesse.m
%This routine evaluates the Hessian of the Rosenbrock function

function H = hesse(x)
H(1,1) = 1200*x(1)^2-400*x(2)+2;
H(1,2) = -400*x(1);
H(2,1) = -400*x(1);
H(2,2) = 200;
end