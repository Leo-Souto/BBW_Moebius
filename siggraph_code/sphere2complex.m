function P_complex = sphere2complex(P_sphere)
P_stereo = [2*P_sphere(:,2)./(P_sphere(:,1)+1) 2*P_sphere(:,3)./(P_sphere(:,1)+1)];
P_complex = P_stereo(:,1) + P_stereo(:,2)*1i;
end