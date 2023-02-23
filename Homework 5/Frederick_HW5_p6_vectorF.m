% % This code integrates on [x0,xf] a model equation 
%Vector function

function z = Frederick_HW5_p6_vectorF(t, y)
    z(1) = y(2); 
    z(2) = -y(1);
    
    z = transpose(z);

end 