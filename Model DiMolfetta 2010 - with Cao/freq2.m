clear
clc

Rc = 0.0398 - 0.005;
Ls = 0.001025;
Rs = 0.8738;
Cs = 2.896;
Cr = 4.0;
Cao = 0.2;
Rc_line = 0.005;

z1 = Rs*Ls*Cr*Cs;
z2 = Rc*Rs*Cr*Cs + Ls*Cr + Ls*Cs;
z3 = Cr*Rs + Rc*Cr + Rc*Cs;
z4 = 1;
p1 = Rs*Cr*Cs;
p2 = Cr + Cs;

r1 = z1*Cao;
r2 = z2*Cao;
r3 = p1 + z3*Cao;
r4 = p2 + Cao;
l1 = Rc_line*r1;
l2 = z1 + Rc_line*r2;
l3 = z2 + Rc_line*r3;
l4 = z3 + Rc_line*r4;

s = tf([1 0],1);
gs1 = (z1*s^3 + z2*s^2 + z3*s + 1)/(p1*s^2 + p2*s);
% bode(gs1)
% hold on
gs2 = (l1*s^4 + l2*s^3 + l3*s^2 + l4*s + 1)/(r1*s^4 + r2*s^3 + r3*s^2 + r4*s);
bode(gs2)
% axis([0.1 60 -30 10])

