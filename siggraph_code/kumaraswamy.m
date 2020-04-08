function result = kumaraswamy(x,a,b,z_max, z_min)
if nargin <= 3
   z_min = 0;
   z_max = 1;
end
y = (x - z_min)./(z_max - z_min);
y = max(min(1,y),0);
result = a.*b.*y.^(a-1).*(1-y.^a).^(b-1);
end