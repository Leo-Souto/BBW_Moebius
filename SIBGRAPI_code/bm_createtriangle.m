function [xp, yp, zp] = bm_createtriangle(P1, P2, flag)

[n1,n2,n3] = bm_planecoeffs(P1,P2);
n = [n1 n2 n3];
[p3x,p3y,p3z] = bm_planecoeffs(P1,n);
[p4x,p4y,p4z] = bm_planecoeffs(P2,n);
P3 = [p3x p3y p3z];
P4 = [p4x p4y p4z];
if flag == 1
    if P4(2) ~= 0 && P4(1) ~= 0 
        if (P3(1)*P4(2)-P3(2)*P4(1)) ~= 0
            t_intercept = ((P1(2)-P2(2))/P4(2) - (P1(1)-P2(1))/P4(1));
            t_intercept = t_intercept*(P4(1)*P4(2)/(P3(1)*P4(2)-P3(2)*P4(1)));
            xp = P1(1) + t_intercept*P3(1);
            yp = P1(2) + t_intercept*P3(2);
            zp = P1(3) + t_intercept*P3(3);
        else
            [xp, yp, zp] = bm_createtriangle(P1, P2, 2);
        end
    else
        [xp, yp, zp] = bm_createtriangle(P1, P2, 2);
    end
elseif flag == 2
    if P4(3) ~= 0 && P4(1) ~= 0 
        if (P3(1)*P4(3)-P3(3)*P4(1)) ~= 0
            t_intercept = ((P1(3)-P2(3))/P4(3) - (P1(1)-P2(1))/P4(1));
            t_intercept = t_intercept*(P4(1)*P4(3)/(P3(1)*P4(3)-P3(3)*P4(1)));
            xp = P1(1) + t_intercept*P3(1);
            yp = P1(2) + t_intercept*P3(2);
            zp = P1(3) + t_intercept*P3(3);
        else
            [xp, yp, zp] = bm_createtriangle(P1, P2, 3);
        end
    else
        [xp, yp, zp] = bm_createtriangle(P1, P2, 3);
    end
elseif flag == 3
    if P4(3) ~= 0 && P4(2) ~= 0 
        if (P3(2)*P4(3)-P3(3)*P4(2)) ~= 0
            t_intercept = ((P1(3)-P2(3))/P4(3) - (P1(2)-P2(2))/P4(2));
            t_intercept = t_intercept*(P4(2)*P4(3)/(P3(2)*P4(3)-P3(3)*P4(2)));
            xp = P1(1) + t_intercept*P3(1);
            yp = P1(2) + t_intercept*P3(2);
            zp = P1(3) + t_intercept*P3(3);
        else
            [xp, yp, zp] = bm_createtriangle(P1, P2, 4);
        end
    else
        [xp, yp, zp] = bm_createtriangle(P1, P2, 4);
    end
else
    error('problema mal definido')
end

% t = linspace(0,2*pi,300);
% x = cos(t)*P1(1) + sin(t)*P2(1);
% y = cos(t)*P1(2) + sin(t)*P2(2);
% z = cos(t)*P1(3) + sin(t)*P2(3);
% 
% plot3(x,y,z)

% grid on
% hold on
% sphere
% plot3([P1(1) P2(1) xp P1(1)],[P1(2) P2(2) yp P1(2)],[P1(3) P2(3) zp P1(3)])
% plot3([P1(1) P2(1) xp],[P1(2) P2(2) yp],[P1(3) P2(3) zp],'*')
% hold off
% axis equal
% cameratoolbar

end