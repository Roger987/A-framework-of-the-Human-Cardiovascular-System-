function y = pumpSpeed(passo, tf, tc, min, max)

    t = 0:passo:tf;
    N = tf/passo;

    y = [];
    for i=1:N+1
        tn = rem(t(i),tc);
        if tn <= tc/2
            y = [y (max+0.005*i)*sin(2*pi*(1/tc)*tn)];
%             y = [y (max)*sin(2*pi*(1/tc)*tn)];
        else
            y = [y min*sin(2*pi*(1/tc)*tn)];
%             y = [y (min+0.01*i)*sin(2*pi*(1/tc)*tn)];
        end

    end

end