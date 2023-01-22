function [Plved_target, Qlvad_target] = target(x1,y1,CL,Plved,SF)
    
    if SF ~= 1
        ControlLine = CL.*SF;
    else
       ControlLine = CL; 
    end
    slope = -1.96;
    y = slope.*(Plved - x1) + y1;
    epsilon = 1e-6;
    idx = find(y - ControlLine < epsilon, 1); %// Index of coordinate in array
    Plved_target = Plved(idx);
    Qlvad_target = y(idx);
    
end