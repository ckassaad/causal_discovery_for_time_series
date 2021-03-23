function showMatrix(mat)
% figure
mat = (mat - min(mat(:))) ./ (max(mat(:)) - min(mat(:)));
imshow(0.9*mat, 'InitialMagnification', 800);