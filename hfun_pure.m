function F = hfun_pure(T,h)
A_pure = [2.502005e-1,-5.1536879e-5,6.5519486e-8,-6.7178376e-12,-1.5128259e-14,7.6215767e-18,-1.4526770e-21,1.0115540e-25];
href_pure = -1.7558886;
F = href_pure + A_pure(1)*T+A_pure(2)/2*T.^2+A_pure(3)/3*T.^3+A_pure(4)/4*T.^4+A_pure(5)/5*T.^5+A_pure(6)/6*T.^6+A_pure(7)/7*T.^7+A_pure(8)/8*T.^8-h;
end