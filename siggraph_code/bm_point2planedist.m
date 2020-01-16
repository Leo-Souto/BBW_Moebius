function [distance] = bm_point2planedist(P,a, b, c, d)
distance = abs(a.*P(:,1) + b.*P(:,2) + c.*P(:,3) + d)./sqrt(a.*a + b.*b +c.*c);
end