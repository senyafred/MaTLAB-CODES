% This is the discontinuous function for Problem 3:

function f=frederick_fun4_hw2_p3(t,y)
global g k1 k2

if t>=0 && t<=2
f= -g - k1*y;
else
f = -g - k2*y;
end
