function [elv, erv, ea] = elastance(passo,tf,tc)

    TPB = 0.92*tc; % P-wave beginning time
    TPE = 1.01*tc; % P-wave ending time

    TT = 0.30*tc; % T-wave peak time
%     TT = 0.35*tc; % T-wave peak time
    TTE = 0.45*tc; % T-wave ending time

    aa = []; % Atrial activation function
    av = []; % Ventricular activation function

    for t = 0:passo:tf

        aa_ = 0;
        av_ = 0;
        tn = rem(t,tc);

        % Atrium
        if tn >= TPB && tn < TPE
            aa_ = 1 - cos(((tn - TPB)/(TPE - TPB))*2*pi);
        end

        % Ventricle
        if tn >=0 && tn < TT
            av_ = 1 - cos((tn/TT)*pi);
        elseif tn >= TT && tn < TTE
            av_ = 1 + cos(((tn - TT)/(TTE - TT))*pi);
        end

        aa = [aa aa_];
        av = [av av_];

    end

    % LEFT AND RIGHT ATRIUM
    e_min = 0.15;
    e_max = 0.25;
    ea = e_min + ((e_max-e_min)/2)*aa;

    % LEFT VENTRICLE
    elv_min = 0.1;
    elv_max = 2.5;
    elv = elv_min + ((elv_max - elv_min)/2)*av;

    % RIGHT VENTRICLE
%     erv_min = 0.025;
    erv_min = 0.0225;
    erv_max = 1.15;
    erv = erv_min + ((erv_max - erv_min)/2)*av;

%     PLOT
%     t = 0:passo:tf;
%     plot(t,elv,t,ea,'-.',t,erv,'--');
%     legend('E(lv)', 'E(la)', 'E(rv)');

end
