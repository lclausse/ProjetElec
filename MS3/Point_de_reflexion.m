function [I,theta, distanceReflexion, check]=Point_de_reflexion(n,V0,P0,P1)

% Inputs: 
%       n: normal vector of the Plane 
%       V0: any point that belongs to the Plane 
%       P0: position of the transmitter
%       P1: position of the receiver
%
%Outputs:
%      I    is the reflexion  point 
%     Check is an indicator:
%      0 => disjoint (no intersection)
%      1 => the plane intersects P0P1 in the unique point I
%      2 => the segment lies in the plane
%      3=>the intersection lies outside the segment P0P1
%
%
%p = n(1)*x+n(2)*y+n(3)*z+(-V0(1) -V0(2) -V0(3))=0
%s<=>   x = n(1)*t + P1(1)
%       y = n(2)*t + P1(2)
%       z = n(3)*t + P1(3)
%  
%M= (P1(1)+n(1)*t, P1(2) + n(2)*t, P1(3) + n(3)*t)
%P2 = (2 * x_M - x_P1, 2 * y_M - y_P1, 2 * z_M - z_P1)



d = -V0(1)-V0(2)-V0(3);
syms t
T = solve( (n(1)*(n(1)* t +P1(1)) + n(2)*(n(2)* t +P1(2)) + n(3)*(n(3)* t +P1(3)) + d == 0), t );

PM =[P1(1)+n(1)*T, P1(2) + n(2)*T, P1(3) + n(3)*T];
P2 = [2*PM(1)-P1(1),2*PM(2)-P1(2),2*PM(3)-P1(3)];
I=[0 0 0];
u = P2-P0;
w = P0 - V0;
D = dot(n,u);
N = -dot(n,w);
check=0;
if abs(D) < 10^-7        % The segment is parallel to plane
        if N == 0           % The segment lies in plane
            check=2;
            return
        else
            check=0;       %no intersection
            return
        end
end
%compute the intersection parameter
sI = N / D;
I = P0+ sI.*u;
if (sI < 0 || sI > 1)
    check= 3;          %The intersection point  lies outside the segment, so there is no intersection
else
    check=1;
end
v = P0-I;

theta = atan2d(norm(cross(n,v)),dot(n,v));
distanceReflexion = norm(P0-I)+norm(P1-I);
[I,theta, distanceReflexion, check];
end
