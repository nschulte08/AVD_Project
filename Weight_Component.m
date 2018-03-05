function [W_nacelle] = Weight_Component( K_ng, N_Lt, N_w, N_z, W_ec, N_en, S_n )
% This script determines the components of the vehicle empty weight 

W_nacelle = 0.6724 * K_ng * N_Lt^0.1 * N_w^0.294 * N_z^0.119 * W_ec^0.611 * N_en^0.984 * S_n^0.224;








end