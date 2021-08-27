function xplus = g_syst(x,lambda)


% state
x1 = x(1);
x2 = x(2);

xplus = [-x1 ; -lambda*x2];
end