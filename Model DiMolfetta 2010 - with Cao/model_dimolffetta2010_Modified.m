clear
clc

passo = 0.0001;
tf = 10;

HR = 72;
tc = 60/HR;

t = 0:passo:tf;
N = tf/passo;

% PARAMETERS

% Sistemic arterial section
Rcs = 0.0495;
% Ls = 0.0000736;
Ls = 0.0000736;
Cas = 3.47;
% Ras = 0.7905;
Ras = 2;
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
Vlv0 = 5.0;
Vla0 = 4.0;
Vrv0 = 10.0;
% Vrv0 = 20.0;
Vra0 = 4.0;

% DIODES
Dm = 0;
Da = 0;
Dt = 0;
Dp = 0;

Dm_v = [Dm];
Da_v = [Da];
Dt_v = [Dt];
Dp_v = [Dp];

% ELASTANCES 
[Elv, Erv, Ea] = elastance(passo,tf,tc);
% [Elv, Erv] = elastance_simaan(passo,tf,tc);

% VARIABLES

Vla = zeros(1,N);
Vlv = zeros(1,N);
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

LEDV = [];
LESV = [];
REDV = [];
RESV = [];

% INITIAL CONDITIONS

% Vla(1) = 48;
% Vlv(1) = 140; 
% Qa(1) = 0; 
% Pas(1) = 75; 
% Pvs(1) = 9.25; 
% Vra(1) = 17; 
% Vrv(1) = 120; 
% Pap(1) = 9; 
% Qp(1) = 0; 
% Pvp(1) = 7.4;

Vla(1) = 68;
Vlv(1) = 125;
Pao(1) = 80;
Qa(1) = 0; 
Pas(1) = 72; 
Pvs(1) = 8.8; 
Vra(1) = 20; 
Vrv(1) = 125; 
Pap(1) = 12; 
Qp(1) = 0; 
Pvp(1) = 10;

Pla(1) = Ea(1)*(Vla(1) - Vla0);
Plv(1) = Elv(1)*(Vlv(1) - Vlv0);
Pra(1) = Ea(1)*(Vra(1) - Vra0);
Prv(1) = Erv(1)*(Vrv(1) - Vrv0);

x = [Vla(1), Vlv(1), Pao(1), Qa(1), Pas(1), Pvs(1), Vra(1), Vrv(1), Pap(1), Qp(1), Pvp(1)]';

for i=1:N
    
    tn = rem(t(i),tc);
    
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
        elseif Da_v(i) == 0 && Da_v(i-1) == 1
            LESV = [LESV Vlv(i)];
            RESV = [RESV Vrv(i)];
        end
    end
    
    % EQUATIONS
    
    a11 = -((1/Rvp)+(Dm/Rm))*Ea(i);
    a12 = (Dm/Rm)*Elv(i);
    a110 = 1/Rvp;
    
    a21 = (Dm/Rm)*Ea(i);
    a22 = -((Dm/Rm)+(Da/Ra))*Elv(i);
    a23 = Da/Ra;
    
    a32 = (Da/Ra)*Elv(i);
    a33 = -Da/Ra;
    
    a44 = -Rcs;
    
    a55 = -1/Ras;
    a56 = 1/Ras;
    
    a65 = 1/Ras;
    a66 = -((1/Ras)+(1/Rvs));
    a67 = Ea(i)/Rvs;
    
    a76 = 1/Rvs;
    a77 = -((1/Rvs)+(Dt/Rt))*Ea(i);
    a78 = (Dt/Rt)*Erv(i);
    
    a87 = (Dt/Rt)*Ea(i);
    a88 = -(Dt/Rt)*Erv(i);
    
    a99 = -1/Rap;
    a911 = 1/Rap;
    
    a108 = Erv(i);
    a1010 = -(Rp + Rcp);
    
    a111 = Ea(i)/Rvp;
    a119 = 1/Rap;
    a1111 = -((1/Rap)+(1/Rvp));
    
    A = zeros(11,11);
    % A =     [Vla   Vlv  Pao    Qa  Pas  Pvs  Vra   Vrv    Pap     Qp   Pvp]
    A(1,:)  = [a11 , a12,   0 ,   0,   0,   0,   0,    0,    0,     0,  a110];         % Vla  
    A(2,:)  = [a21 , a22,  a23,   0,   0,   0,   0,    0,    0,     0,     0];         % Vlv
    A(3,:)  = [ 0  , a32,  a33,  -1,   0,   0,   0,    0,    0,     0,     0]/Cao;     % Pao
    A(4,:)  = [ 0  ,  0 ,   1 , a44,  -1,   0,   0,    0,    0,     0,     0]/Ls;      % Qa
    A(5,:)  = [ 0  ,  0 ,   0 ,   1, a55, a56,   0,    0,    0,     0,     0]/Cas;     % Pas
    A(6,:)  = [ 0  ,  0 ,   0,    0, a65, a66, a67,    0,    0,     0,     0]/Cvs;     % Pvs
    A(7,:)  = [ 0  ,  0 ,   0,    0,   0, a76, a77,  a78,    0,     0,     0];         % Vra
    A(8,:)  = [ 0  ,  0 ,   0,    0,   0,   0, a87,  a88,    0,    -1,     0];         % Vrv
    A(9,:)  = [ 0  ,  0 ,   0,    0,   0,   0,   0,    0,  a99,     1,  a911]/Cap;     % Pap
    A(10,:) = [ 0  ,  0 ,   0,    0,   0,   0,   0, a108,   -1, a1010,     0]*(Dp/Lp); % Qp
    A(11,:) = [a111,  0 ,   0,    0,   0,   0,   0,    0, a119,     0, a1111]/Cvp;     % Pvp
 

    B = [              ((1/Rvp)+(Dm/Rm))*(Ea(i)*Vla0) - ((Dm/Rm)*Elv(i)*Vlv0);   % Vlv
                      -((Dm/Rm)*(Ea(i)*Vla0)) + ((Dm/Rm)+(Da/Ra))*Elv(i)*Vlv0;   % Vla
                                                   -(Da/(Ra*Cao))*Elv(i)*Vlv0;   % Pao
                                                                            0;   % Qa
                                                                            0;   % Pas
                                                  -(1/(Rvs*Cvs))*(Ea(i)*Vra0);   % Pvs
                       ((1/Rvs)+(Dt/Rt))*(Ea(i)*Vra0) - (Dt/Rt)*(Erv(i)*Vrv0);   % Vra
                                         (Dt/Rt)*(- Ea(i)*Vra0 + Erv(i)*Vrv0);   % Vrv 
                                                                            0;   % Pap
                                                        -(1/Lp)*(Erv(i)*Vrv0);   % Qp
                                                  -(1/(Rvp*Cvp))*(Ea(i)*Vla0)];  % Pvp
    
    if Dp == 0
        x(10) = 0;
        B(10) = 0;
    end

    x = runkut4(passo, A, x, B);
    
    Vla(i+1) = x(1);
    Vlv(i+1) = x(2);
    Pao(i+1) = x(3);
    Qa(i+1)  = x(4);
    Pas(i+1) = x(5);
    Pvs(i+1) = x(6);
    Vra(i+1) = x(7);
    Vrv(i+1) = x(8);
    Pap(i+1) = x(9);
    Qp(i+1)  = x(10);    
    Pvp(i+1) = x(11);
    
    Pla(i+1) = Ea(i+1)*(Vla(i+1) - Vla0);
    Plv(i+1) = Elv(i+1)*(Vlv(i+1) - Vlv0);
    Pra(i+1) = Ea(i+1)*(Vra(i+1) - Vra0);
    Prv(i+1) = Erv(i+1)*(Vrv(i+1) - Vrv0);
       
end

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

