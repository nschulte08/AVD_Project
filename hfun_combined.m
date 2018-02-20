function F = hfun_combined(T,h_f,f)
A_pure = [2.502005e-1,-5.1536879e-5,6.5519486e-8,-6.7178376e-12,-1.5128259e-14,7.6215767e-18,-1.4526770e-21,1.0115540e-25];
A_vit = [7.3816638e-2,1.2258630e-3,-1.3771901e-6,9.9686793e-10,-4.2051104e-13,1.0212913e-16,-1.3335668e-20,7.2678710e-25];
href_pure = -1.7558886;%[Btu/lbm]
href_vit = 30.58153;%[Btu/(lbm ? °R)]
%pg 91 "Elements of Propulsion"%
F = -h_f+href_pure + A_pure(1)*T+A_pure(2)/2*T.^2+A_pure(3)/3*T.^3+A_pure(4)/4*T.^4+A_pure(5)/5*T.^5+A_pure(6)/6*T.^6+A_pure(7)/7*T.^7+A_pure(8)/8*T.^8+f*(href_vit+A_vit(1)*T+A_vit(2)/2*T.^2+A_vit(3)/3*T.^3+A_vit(4)/4*T.^4+A_vit(5)/5*T.^5+A_vit(6)/6*T.^6+A_vit(7)/7*T.^7+A_vit(8)/8*T.^8);%pg 91 "Elements of Propulsion"
end