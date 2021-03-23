%Dino Sejdinovic, 2013
%D. Sejdinovic, A. Gretton and W. Bergsma.  A KERNEL TEST FOR THREE-VARIABLE INTERACTIONS, 2013.

%---median heuristic bandwidth selection

function [sig] = median_heur(Z)

size1=size(Z,1);
    if size1>100 %choose 100 random samples from Z if it contains
        %more than 100 points
        [~,ind]=sort(rand(1,size1));
        Zmed = Z(ind(1:100),:);
        size1 = 100;
    else
      Zmed = Z;
    end
    G = sum((Zmed.*Zmed),2);
    Q = repmat(G,1,size1);
    R = repmat(G',size1,1);
    dists = Q + R - 2*(Zmed*Zmed');
    dists = dists-tril(dists);
    dists = reshape(dists,size1^2,1);
    sig = sqrt(0.5*median(dists(dists>0)));
end

