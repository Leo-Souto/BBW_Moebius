function [a, b, c] = bm_planecoeffs(P1,P2)

a = P1(:,2).*P2(:,3) - P1(:,3).*P2(:,2);
b = P1(:,3).*P2(:,1) - P1(:,1).*P2(:,3);
c = P1(:,1).*P2(:,2) - P1(:,2).*P2(:,1);
end