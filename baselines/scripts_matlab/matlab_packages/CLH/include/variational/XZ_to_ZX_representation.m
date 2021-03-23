%% XZ_to_ZX_representation
%
% Summary: Inverse of ZX_to_XZ_representation


function A_ZX = XZ_to_ZX_representation(A_XZ, K_X)

K = size(A_XZ, 1);
K_Z = K - K_X;

A_ZX = A_XZ( [ (K_X + 1):K , 1:K_X ] , [ (K_X + 1):K , 1:K_X ] );
