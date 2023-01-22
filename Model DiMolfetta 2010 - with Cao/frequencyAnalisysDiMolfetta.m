% clear
% clc

% Without Cao

Rcs = 0.0495;
Ls = 0.0000736;
Ras = 0.7905;
Cas = 3.47;
Cvs = 84;

% Cao = 0.08;

z1 = Ras*Ls*Cvs*Cas;
z2 = Rcs*Ras*Cvs*Cas + Ls*Cvs + Ls*Cas;
z3 = Ras*Cvs + Rcs*Cvs + Rcs*Cas;
p1 = Ras*Cvs*Cas;
p2 = Cvs + Cas;

% r1 = z1*Cao;
% r2 = z2*Cao;
% r3 = p1 + z3*Cao;
% r4 = p2 + Cao;
% l1 = Rcs_line*r1;
% l2 = z1 + Rcs_line*r2;
% l3 = z2 + Rcs_line*r3;
% l4 = z3 + Rcs_line*r4;

% wf = 10*2*pi;
w = 0.4:0.001:10*2*pi;
Z0_jw = (-z1*j.*w.^3 - z2.*w.^2 + z3*j.*w + 1)./(-p1.*w.^2 + p2*j.*w);
subplot(2,1,1)
hold on
plot(w,abs(Z0_jw))
ylabel('Magnitude ()')
subplot(2,1,2)
hold on
plot(w,phase(Z0_jw))
ylabel('Phase (radians)')
xlabel('Frequency (rad/s)')

% hold on
% sys = tf([z1,z2,z3,1],[p1,p2,0])
% bode(sys)
% hold on

% Zm_jw = (l1.*w.^4 - l2*j.*w.^3 - l3.*w.^2 + l4*j.*w + 1)./(r1.*w.^4 - r2*j.*w.^3 - r3.*w.^2 + r4*j.*w);
% plot(w,abs(Z0_jw),w,abs(Zm_jw))
% plot(w,abs(Zm_jw))

%%

% With Cao and Rc'

clear
clc

Rcs_line = 0.005;
Rcs = 0.0495 - Rcs_line;
Ls = 0.0000736;
Ras = 0.7905;
Cas = 3.47;
Cvs = 84;

Cao = 0.08;

z1 = Ras*Ls*Cvs*Cas;
z2 = Rcs*Ras*Cvs*Cas + Ls*Cvs + Ls*Cas;
z3 = Ras*Cvs + Rcs*Cvs + Rcs*Cas;
p1 = Ras*Cvs*Cas;
p2 = Cvs + Cas;

r1 = z1*Cao;
r2 = z2*Cao;
r3 = p1 + z3*Cao;
r4 = p2 + Cao;
l1 = Rcs_line*r1;
l2 = z1 + Rcs_line*r2;
l3 = z2 + Rcs_line*r3;
l4 = z3 + Rcs_line*r4;

% wf = 10*2*pi;
w = 0.4:0.001:60;
% Z0_jw = (-z1*j.*w.^3 - z2.*w.^2 + z3*j.*w + 1)./(-p1.*w.^2 + p2*j.*w);
% plot(w,abs(Z0_jw))

hold on
Zm_jw = (l1.*w.^4 - l2*j.*w.^3 - l3.*w.^2 + l4*j.*w + 1)./(r1.*w.^4 - r2*j.*w.^3 - r3.*w.^2 + r4*j.*w);
% plot(w,abs(Z0_jw),w,abs(Zm_jw))
plot(w,abs(Zm_jw))
hold on

%% 

% With Cao

clear
clc

Rcs_line = 0.0;
Rcs = 0.0495 - Rcs_line;
Ls = 0.0000736;
Ras = 0.7905;
Cas = 3.47;
Cvs = 84;

Cao = 0.08;

z1 = Ras*Ls*Cvs*Cas;
z2 = Rcs*Ras*Cvs*Cas + Ls*Cvs + Ls*Cas;
z3 = Ras*Cvs + Rcs*Cvs + Rcs*Cas;
p1 = Ras*Cvs*Cas;
p2 = Cvs + Cas;

r1 = z1*Cao;
r2 = z2*Cao;
r3 = p1 + z3*Cao;
r4 = p2 + Cao;
l1 = Rcs_line*r1;
l2 = z1 + Rcs_line*r2;
l3 = z2 + Rcs_line*r3;
l4 = z3 + Rcs_line*r4;

% wf = 10*2*pi;
w = 0.4:0.001:10*2*pi;
% Z0_jw = (-z1*j.*w.^3 - z2.*w.^2 + z3*j.*w + 1)./(-p1.*w.^2 + p2*j.*w);
% plot(w,abs(Z0_jw))
hold on
Zm_jw = (l1.*w.^4 - l2*j.*w.^3 - l3.*w.^2 + l4*j.*w + 1)./(r1.*w.^4 - r2*j.*w.^3 - r3.*w.^2 + r4*j.*w);
% plot(w,abs(Z0_jw),w,abs(Zm_jw))
subplot(2,1,1)
plot(w,abs(Zm_jw))

hold on
subplot(2,1,2)
plot(w,phase(Zm_jw))

% hold on
% sys = tf([l1,l2,l3,l4,1],[r1,r2,r3,r4,0])
% bode(sys)
