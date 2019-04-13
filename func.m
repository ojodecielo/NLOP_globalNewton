%file name:  func.m
%This routine evaluates the Rosenbrock function.

function y = func(x)
y = 100*(x(2) - x(1)^2)^2 + (1-x(1))^2;
end