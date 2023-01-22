function [elv, erv] = elastance_simaan(passo,tf,tc)

    t = 0:passo:tf;
    N = tf/passo;
    t_max = 0.2 + 0.15*tc;

    for i = 1:N+1
        tn(i) = mod(i*passo,tc)/t_max;
%         tn(i) = rem(t(i),tc);
        t1 = (tn(i)/0.7)^1.9;
        t2 = 1 + t1;
        t3 = 1 + (tn(i)/1.17)^21.9;
        En(i) = 1.55*(t1/t2)*(1/t3);
    end   
    
    Elv_max = 2.2;
    Elv_min = 0.05;
    Erv_max = 1.15;
    Erv_min = 0.025;
    
    elv = (Elv_max - Elv_min)*En + Elv_min;
    erv = (Erv_max - Erv_min)*En + Erv_min;

%     PLOT
%     t = 0:passo:tf;
%     plot(t,elv,t,ea,'-.',t,erv,'--');
%     legend('E(lv)', 'E(la)', 'E(rv)');

end