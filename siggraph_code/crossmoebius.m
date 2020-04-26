function [a,b,c,d] = crossmoebius(p1, p2, p3, p1_n, p2_n, p3_n)

a = det([p1*p1_n p1_n 1;
         p2*p2_n p2_n 1;
         p3*p3_n p3_n 1]);     
b = det([p1*p1_n p1 p1_n;
         p2*p2_n p2 p2_n;
         p3*p3_n p3 p3_n]);
c = det([p1 p1_n 1;
         p2 p2_n 1;
         p3 p3_n 1]);
d = det([p1*p1_n p1 1;
         p2*p2_n p2 1;
         p3*p3_n p3 1]);

end