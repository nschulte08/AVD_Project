function F = phifun_combined(T,phi,f)
A_pure = [2.502005e-1,-5.1536879e-5,6.5519486e-8,-6.7178376e-12,-1.5128259e-14,7.6215767e-18,-1.4526770e-21,1.0115540e-25];
A_vit = [7.3816638e-2,1.2258630e-3,-1.3771901e-6,9.9686793e-10,-4.2051104e-13,1.0212913e-16,-1.3335668e-20,7.2678710e-25];
phi_ref_pure = 0.0454323;
phi_ref_vit = 0.6483398;
F = -phi+phi_ref_pure + A_pure(1)*log(T)+A_pure(2)*T+A_pure(3)/2*T.^2+A_pure(4)/3*T.^3+A_pure(5)/4*T.^4+A_pure(6)/5*T.^5+A_pure(7)/6*T.^6+A_pure(8)/7*T.^7+f*(phi_ref_vit + A_vit(1)*log(T)+A_vit(2)*T+A_vit(3)/2*T.^2+A_vit(4)/3*T.^3+A_vit(5)/4*T.^4+A_vit(6)/5*T.^5+A_vit(7)/6*T.^6+A_vit(8)/7*T.^7);%eq. 2.640 Elements of Propulsion?
end