function F = phifun_vit(T,phi)
A_vit = [7.3816638e-2,1.2258630e-3,-1.3771901e-6,9.9686793e-10,-4.2051104e-13,1.0212913e-16,-1.3335668e-20,7.2678710e-25];
phi_ref_vit = 0.6483398;
F = -phi+phi_ref_vit + A_vit(1)*ln(T)+A_vit(2)*T+A_vit(3)/2*T.^2+A_vit(4)/3*T.^3+A_vit(5)/4*T.^4+A_vit(6)/5*T.^5+A_vit(7)/6*T.^6+A_vit(8)/7*T.^7;
end