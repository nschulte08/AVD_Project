function F = phifun_pure(T,phi)
A_pure = [2.502005e-1,-5.1536879e-5,6.5519486e-8,-6.7178376e-12,-1.5128259e-14,7.6215767e-18,-1.4526770e-21,1.0115540e-25];
phi_ref_pure = 0.0454323;
F = -phi+phi_ref_pure + A_pure(1)*ln(T)+A_pure(2)*T+A_pure(3)/2*T.^2+A_pure(4)/3*T.^3+A_pure(5)/4*T.^4+A_pure(6)/5*T.^5+A_pure(7)/6*T.^6+A_pure(8)/7*T.^7;
end