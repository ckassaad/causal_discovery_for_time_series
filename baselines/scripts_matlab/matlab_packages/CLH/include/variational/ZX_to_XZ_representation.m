%% Transfrom Mingming's representation of A, i.e. according to order (Z,X), to our representation, i.e. according to order (X,Z)
% Input:
% * A_ZX: the matrix in (Z,X)-from
% * K_X: dimension of X, i.e. number of observed components

function A_XZ = ZX_to_XZ_representation(A_ZX, K_X)

K = size(A_ZX, 1);
K_Z = K - K_X;

A_XZ = A_ZX( [ (K_Z + 1):K , 1:K_Z ] , [ (K_Z + 1):K , 1:K_Z ] );
