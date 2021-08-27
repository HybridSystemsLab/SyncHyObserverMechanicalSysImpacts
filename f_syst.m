function xdot = f_syst(x,gamma,rho)

% state
x1 = x(1);
x2 = x(2);


% differential equations
xdot = [x2 ; gamma-rho*x2];
end