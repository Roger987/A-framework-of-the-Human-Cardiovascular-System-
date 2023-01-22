clear
clc

Rc = 0.0398 - 0.005;
Ls = 0.001025;
Rs = 0.8738;
Cs = 2.896;
Cr = 4.0;

Cao = 0.2;
Rc_line = 0.005;
Rc1 = Rc - Rc_line;


z1 = Rs*Ls*Cr*Cs;
z2 = Rc*Rs*Cr*Cs + Ls*Cr + Ls*Cs;
z3 = Rs*Cr + Rc*Cr + Rc*Cs;
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

% wf = 10*2*pi;
w = 0.4:0.001:60;
Z0_jw = (-z1*j.*w.^3 - z2.*w.^2 + z3*j.*w + 1)./(-p1.*w.^2 + p2*j.*w);
% plot(w,abs(Z0_jw))

Zm_jw = (l1.*w.^4 - l2*j.*w.^3 - l3.*w.^2 + l4*j.*w + 1)./(r1.*w.^4 - r2*j.*w.^3 - r3.*w.^2 + r4*j.*w);
% plot(w,abs(Z0_jw),w,abs(Zm_jw))
plot(w,abs(Zm_jw))

% w = logspace(0,2);
% 
% h = freqs(num1,den1,w);
% mag = abs(h);
% phase = angle(h);
% 
% subplot(2,1,1)
% loglog(w,mag)
% grid on
% xlabel('Frequency (rad/s)')
% ylabel('Magnitude')
% 
% subplot(2,1,2)
% semilogx(w,phase)
% grid on
% xlabel('Frequency (rad/s)')
% ylabel('Phase')


% zm = tf(num1,den1)
% 
% % num1 = [R3*L*C1*C2, R2*R3*C1*C2 + L*C1 + L*C2, R2*C1 + R2*C2 + R3*C2, 1];
% % den1 = [0, R3*C1*C2, C1 + C2, 0];
% % zm = tf(num1,den1)
% 
% % num2 = [Ls, Ra+Rcs];
% % den2 = 1;
% % zo = tf(num2,den2);
% 
% % figure()
% % stepplot(zm)
% 
% figure()
% bode(zm)