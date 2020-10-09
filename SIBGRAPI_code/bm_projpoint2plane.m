function [P_projetado] = bm_projpoint2plane(P,a,b,c,d)

% dist = (P(:,1)*a+P(:,2)*b+P(:,3)*c) + d;
% P_projetado = [P(:,1)-dist*a P(:,2)-dist*b P(:,3)-dist*c];

D = 0; E = 0; F = -d/c;
t = (a*D-a*P(:,1) + b*E-b*P(:,2) + c*F-c*P(:,3))/(a^2+b^2+c^2);
P_projetado = P+t*[a b c];

end