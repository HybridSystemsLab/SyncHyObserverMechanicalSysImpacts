function [value] = C_syst(x) 


x1 = x(1);

if x1 >= 0
    value = 1;
else
    value = 0;
end
end