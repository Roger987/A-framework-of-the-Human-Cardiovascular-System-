clear
clc

passo = 0.0001;
tf = 5;
t = passo:passo:tf;
N = tf/passo;

% HR = zeros(1,N);
% HR(1) = 72;

% PARAMETERS

% Sistemic arterial section
Rcs = 0.0495;
Ls  = 0.0000736;
Cas = 3.47;
Ras = 0.7905;
Cao = 0.08;

% Pulmonary arterial section
Rcp = 0.00225;
Lp = 0.000035;
Cap = 4;
Rap = 0.075;

% Sistemic venous section
Cvs = 84;
Rvs = 0.06;

% Pulmonary venous section
Cvp = 5;
Rvp = 0.00075;

% Left heart
Rm = 0.003;
Ra = 0.0075;

% Right heart
Rt = 0.00525;
Rp = 0.00225;

% Constants
Vlv0 = 4.0;
Vla0 = 4.0;
% Vrv0 = 10.0;
Vrv0 = 20.0;
Vra0 = 4.0;

% VAD parameters
Ri = 0.0677;
Ro = 0.0677;
Li = 0.0127;
Lo = 0.0127;
Rk = 0;
alpha = -3.5;
B0 = 0.17070;
B1 = 0.02177;
B2 = -0.000093;
x1 = 1.0;

% DIODES
Dm = 0;
Da = 0;
Dt = 0;
Dp = 0;
Dm_v = Dm;
Da_v = Da;
Dt_v = Dt;
Dp_v = Dp;

% ELASTANCES 
% [Elv, Erv, Ea] = elastance(passo,tf,tc);
Elv = zeros(1,N);
Erv = zeros(1,N);
Ea = zeros(1,N);

% VARIABLES
Vla = zeros(1,N);
Vlv = zeros(1,N);
Qvad = zeros(1,N);
Pao = zeros(1,N);
Qa = zeros(1,N);
Pas = zeros(1,N);
Pvs = zeros(1,N);
Vra = zeros(1,N);
Vrv = zeros(1,N);
Pap = zeros(1,N);
Qp = zeros(1,N);
Pvp = zeros(1,N);

Pla = zeros(1,N);
Plv = zeros(1,N);
Pra = zeros(1,N);
Prv = zeros(1,N);

w = zeros(1,N); % pump speed
Plved = zeros(1,N); % left-ventricular end-diastolic pressure
PlvedTarget = zeros(1,N);
Qvad_mean_v = zeros(1,N);
Qvad_target = zeros(1,N);

Rm_v = zeros(1,N);

LEDV = [];
LESV = [];
REDV = [];
RESV = [];

% INITIAL CONDITIONS
Vla(1) = 68;
Vlv(1) = 125;
Qvad(1) = 0;
Pao(1) = 75;
Qa(1) = 0; 
Pas(1) = 72; 
Pvs(1) = 8.75; 
Vra(1) = 20; 
Vrv(1) = 125; 
Pap(1) = 12; 
Qp(1) = 0; 
Pvp(1) = 10;

Pla(1) = 9.5;
Plv(1) = 13;
Pra(1) = 2;
Prv(1) = 15;

w(1) = 950;

Plved_current = 13;
Plved(1) = Plved_current;

Plved_CL = 0:0.0001:25;
ControlLine = (10.3+(-10.3./(1+(Plved_CL/7).^2.3)));

W_Qvad = Qvad(1);
Qvad_mean = 7.2;
Qvad_mean_v(1) = Qvad_mean;
[PlvedTarget(1),Qvad_target(1)] = target(Plved_current, Qvad_mean, ControlLine, Plved_CL, 1);

Rm_v(1) = 0.003;

% CONTROLLER VARIABLES
error = zeros(1,N);
ierror = zeros(1,N);
Kp = 40;
Ki = 75;

x = [Vla(1), Vlv(1), Qvad(1), Pao(1), Qa(1), Pas(1), Pvs(1), Vra(1), Vrv(1), Pap(1), Qp(1), Pvp(1)]'; 

HR = 90*ones(1,N);

for i=1:N-1
    
    tc = 60/HR(i);
    
    tn = rem(t(i),tc);
    
    [ElvCurrent, ErvCurrent, EaCurrent] = elastancePoint(tn,tc);
    Elv(i+1) = ElvCurrent;
    Erv(i+1) = ErvCurrent;
    Ea(i+1) = EaCurrent;
    
    if Pla(i) >= Plv(i)
        Dm = 1;
    else
        Dm = 0;
    end
    
    if Plv(i) >= Pao(i)
        Da = 1;
    else
        Da = 0;
    end
    
    if Pra(i) >= Prv(i)
        Dt = 1;
    else
        Dt = 0;
    end
    
    if Prv(i) >= Pap(i)
        Dp = 1;
    else
        Dp = 0;
    end
    
    Dm_v = [Dm_v Dm];
    Da_v = [Da_v Da];
    Dt_v = [Dt_v Dt];
    Dp_v = [Dp_v Dp];
    
    if i > 1
        if Dm_v(i) == 0 && Dm_v(i-1) == 1
            LEDV = [LEDV Vlv(i)];
            REDV = [REDV Vrv(i)];
            Plved_current = Plv(i);
            Qvad_mean = mean(W_Qvad);
            W_vad = [];
        elseif Da_v(i) == 0 && Da_v(i-1) == 1
            LESV = [LESV Vlv(i)];
            RESV = [RESV Vrv(i)];
        end
    end

    if Plv(i) <= x1
        Rk = alpha*(Plv(i) - x1);
    else
        Rk = 0;
    end
       
    % EQUATIONS  
    RR = Ri + Ro + Rk + B0;
    LL = Li + Lo + B1;
    
    a11 = -((1/Rvp)+(Dm/Rm))*Ea(i);
    a12 = (Dm/Rm)*Elv(i);
    a112 = 1/Rvp;
    
    a21 = (Dm/Rm)*Ea(i);
    a22 = -((Dm/Rm)+(Da/Ra))*Elv(i);
    a24 = Da/Ra;
    
    a32 = Elv(i);
    a33 = -RR;
    
    a42 = (Da/Ra)*Elv(i);
    a44 = -(Da/Ra);
    
    a55 = -Rcs;
    
    a66 = -1/Ras;
    a67 = 1/Ras;
    
    a76 = 1/Ras;
    a77 = -((1/Ras)+(1/Rvs));
    a78 = Ea(i)/Rvs;
   
    a87 = 1/Rvs;
    a88 = -((1/Rvs)+(Dt/Rt))*Ea(i);
    a89 = (Dt/Rt)*Erv(i);
    
    a98 = (Dt/Rt)*Ea(i);
    a99 = -(Dt/Rt)*Erv(i);
    
    a1010 = -1/Rap;
    a1012 = 1/Rap;
    
    a119 = Erv(i);
    a1111 = -(Rp + Rcp);
    
    a121 = Ea(i)/Rvp;
    a1210 = 1/Rap;
    a1212 = -((1/Rap)+(1/Rvp));
    
    % A =     [Vla   Vlv  Qvad  Pao   Qa  Pas  Pvs    Vra   Vrv    Pap     Qp  Pvp]
    A = zeros(12,12);
    A(1,:)  = [a11 , a12,   0,   0,   0,   0,   0,    0,    0,     0,     0,  a112];         % Vla  
    A(2,:)  = [a21 , a22,  -1, a24,   0,   0,   0,    0,    0,     0,     0,     0];         % Vlv
    A(3,:)  = [ 0  , a32, a33,  -1,   0,   0,   0,    0,    0,     0,     0,     0]/LL;      % Qvad
    A(4,:)  = [ 0  , a42,   1, a44,  -1,   0,   0,    0,    0,     0,     0,     0]/Cao;     % Pao
    A(5,:)  = [ 0  ,  0 ,   0,   1, a55,  -1,   0,    0,    0,     0,     0,     0]/Ls;      % Qa
    A(6,:)  = [ 0  ,  0 ,   0,   0,   1, a66, a67,    0,    0,     0,     0,     0]/Cas;     % Pas
    A(7,:)  = [ 0  ,  0 ,   0,   0,   0, a76, a77,  a78,    0,     0,     0,     0]/Cvs;     % Pvs
    A(8,:)  = [ 0  ,  0 ,   0,   0,   0,   0, a87,  a88,  a89,     0,     0,     0];         % Vra
    A(9,:)  = [ 0  ,  0 ,   0,   0,   0,   0,   0,  a98,  a99,     0,    -1,     0];         % Vrv
    A(10,:) = [ 0  ,  0 ,   0,   0,   0,   0,   0,    0,    0, a1010,     1, a1012]/Cap;     % Pap
    A(11,:) = [ 0  ,  0 ,   0,   0,   0,   0,   0,    0, a119,    -1, a1111,     0]*(Dp/Lp); % Qp
    A(12,:) = [a121,  0 ,   0,   0,   0,   0,   0,   0,    0,  a1210,     0, a1212]/Cvp;     % Pvp
 
    B = [                    ((1/Rvp)+(Dm/Rm))*(Ea(i)*Vla0) - (Dm/Rm)*(Elv(i)*Vlv0); % Vla
                            -(Dm/Rm)*(Ea(i)*Vla0) + ((Dm/Rm)+(Da/Ra))*(Elv(i)*Vlv0); % Vlv
                                              - (B2*(w(i)^2))/LL - (Elv(i)*Vlv0)/LL; % Qvad
                                                         -(Da*Elv(i)*Vlv0)/(Ra*Cao); % Pao
                                                                                  0; % Qa
                                                                                  0; % Pas
                                                            -(Ea(i)*Vra0)/(Rvs*Cvs); % Pvs
                         (((1/Rvs)+(Dt/Rt))*(Ea(i)*Vra0)) - ((Dt/Rt)*(Erv(i)*Vrv0)); % Vra
                                               (Dt/Rt)*(- Ea(i)*Vra0 + Erv(i)*Vrv0); % Vrv
                                                                                  0; % Pap
                                                                  -(Erv(i)*Vrv0)/Lp; % Qp
                                                           -(Ea(i)*Vla0)/(Rvp*Cvp)]; % Pvp
  
    if Dp == 0
        x(11) = 0;
        B(11) = 0;
    end

    x = runkut4(passo, A, x, B);
    
    Vla(i+1)  = x(1);
    Vlv(i+1)  = x(2);
    Qvad(i+1) = x(3);
    Pao(i+1)  = x(4);
    Qa(i+1)   = x(5);
    Pas(i+1)  = x(6);
    Pvs(i+1)  = x(7);
    Vra(i+1)  = x(8);
    Vrv(i+1)  = x(9);
    Pap(i+1)  = x(10);
    Qp(i+1)   = x(11);    
    Pvp(i+1)  = x(12); 
    
    Pla(i+1) = Ea(i+1)*(Vla(i+1) - Vla0);
    Plv(i+1) = Elv(i+1)*(Vlv(i+1) - Vlv0);
    Pra(i+1) = Ea(i+1)*(Vra(i+1) - Vra0);
    Prv(i+1) = Erv(i+1)*(Vrv(i+1) - Vrv0);
    
    % CONTROL
    W_Qvad = [W_Qvad Qvad(i+1)*0.06];
    Qvad_mean_v(i+1) = Qvad_mean;
    Plved(i+1) = Plved_current;
    [PlvedTarget(i+1),Qvad_target(i+1)] = target(Plved_current, Qvad_mean, ControlLine, Plved_CL, 1); % IMMEDIATE STARLING-LIKE
    
    error(i+1)  = Qvad_target(i+1) - Qvad_mean_v(i+1);
    ierror(i+1) = ierror(i) + passo*error(i);
    
    w(i+1) = Kp*error(i+1) + Ki*ierror(i+1); % PI CONTROL
    if i > 120000
        Rm = 0.05;
    end
    Rm_v(i) = Rm;

end
%%
figure
plot(Plved(90000:end),Qvad_mean_v(90000:end),'b','LineWidth',1);
hold on
plot(Plved_CL,ControlLine,'r','LineWidth',1.5);
plot(Plved(99999),Qvad_mean_v(99999),'^-','MarkerFaceColor','green','MarkerEdgeColor','black')
plot(Plved(end),Qvad_mean_v(end),'o-','MarkerFaceColor','green','MarkerEdgeColor','black')
axis([5 10 4.5 7]);
xlabel('P_{lved} (mmHg)','FontSize',13,'FontWeight','bold') 
ylabel('Q_{vad} (L/min)','FontSize',13,'FontWeight','bold') 

%%
% yyaxis left
plot(t,Qvad_mean_v,'b','LineWidth',1.5);
hold on
plot(t,Qvad_target,'-.r','LineWidth',1.5);
leg1 = legend('$\overline{Q_{vad}}$','$\overline{Q_{vad}}^*$');
set(leg1,'Interpreter','latex');
set(leg1,'FontSize',20);
ylabel('Measured Flow and Target Flow (L/min)','FontSize',14,'FontWeight','bold') 
xlabel('Time (s)','FontSize',14,'FontWeight','bold')

%%

% REDV(end)
% RESV(end)
% LEDV(end) 
% LESV(end)

SV_right = (REDV(end) - RESV(end));
SV_left = (LEDV(end) - LESV(end));

CO_right = HR*0.001*SV_right;
CO_left = HR*0.001*SV_left;

ejectionFraction_right = 100*(SV_right/REDV(end));
ejectionFraction_left = 100*(SV_left/LEDV(end));

n_cycles = round(tf/tc);
MAP = zeros(1,n_cycles);
previous_cycle = 1;

for i=1:n_cycles
    MAP(i) = mean(Pas(previous_cycle : i*(N+1)/n_cycles));
    previous_cycle = i*(N+1)/n_cycles;
end

%%
yyaxis left
plot(t,Plved,'b','LineWidth',1.5)
ylabel('Left-Ventricular End-Diastolic Pressure (mmHg)','FontSize',14,'FontWeight','bold') 
xlabel('Time (s)','FontSize',14,'FontWeight','bold')

hold on
yyaxis right
plot(t,Rm_v,'-.r','LineWidth',1.5)

leg1 = legend('$P_{lved}$','$R_{m}$');
ylabel('Mitral Valve Resistance (mmHg s/mL)','FontSize',14,'FontWeight','bold') 

set(leg1,'Interpreter','latex');
set(leg1,'FontSize',20);

ax = gca;
ax.YAxis(1).Color = 'blue';
ax.YAxis(2).Color = 'red';

axis([0 20 0 0.1])

%%
yyaxis left
plot(t,Qvad,'b','LineWidth',1.5)
ylabel('Pump Flow (mL/s)','FontSize',14,'FontWeight','bold') 
xlabel('Time (s)','FontSize',14,'FontWeight','bold')

hold on
yyaxis right
plot(t,w,'r','LineWidth',1.5)
ylabel('Pump Speed (rpm)','FontSize',14,'FontWeight','bold') 
ax = gca;
ax.YAxis(1).Color = 'blue';
ax.YAxis(2).Color = 'red';

leg1 = legend('$Q_{vad}$','$\omega$');
set(leg1,'Interpreter','latex');
set(leg1,'FontSize',20);

