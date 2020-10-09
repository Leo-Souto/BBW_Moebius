function k = bm_intersec_retas(P,P1,P2)
t_1 = (P1(1).*P(:,2) - P1(2).*P(:,1))./((P2(2)-P1(2)).*P(:,1) - (P2(1)-P1(1)).*P(:,2));
if ~iscolumn(t_1)
    t_1 = t_1';
end
t_2 = (P1(1).*P(:,3) - P1(3).*P(:,1))./((P2(3)-P1(3)).*P(:,1) - (P2(1)-P1(1)).*P(:,3));
if ~iscolumn(t_2)
    t_2 = t_2';
end
t_3 = (P1(2).*P(:,3) - P1(3).*P(:,2))./((P2(3)-P1(3)).*P(:,2) - (P2(2)-P1(2)).*P(:,3));
if ~iscolumn(t_3)
    t_3 = t_3';
end

T = horzcat(t_1,t_2,t_3);
T = isfinite(T).*T;
T = max(T,[],2);
T(isnan(T)) = Inf;
k = T;
end