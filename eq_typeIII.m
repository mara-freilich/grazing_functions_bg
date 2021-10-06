function [P,Z] = eq_typeIII(B)
% Function to compute equilibrium phytoplankton and zooplankton
% concentrations for model using type III response in Freilich et al (in 
% review) Biogeosciences
% Input: 
% B: parameters for biological model
N_max = 30;
P = zeros(1,365);
Z = zeros(1,365);
for year_day = 1:365
    x0 = [20 10];
    options1 = odeset('Refine',1,'NonNegative',1);
    t = linspace(0,3650,3650);
    [~,y1]=ode45(@f_typeIII,t,x0,options1);
    
    P(year_day) = y1(end,1);
    Z(year_day) = y1(end,2);
end

    function dydt = f_typeIII(t,y)
        h1 = 20;
        loffset = 270;
        a = 0.5;
        mumax = 0.8;
        N0 = 4;
        [MLD,~,~,~,~,~] = mldmodel(year_day);
        
        yd1 = mod(year_day+loffset,365);
        h_light=20*(0.6*sin(yd1*pi/365*2)+1);
        light=h_light*h1/MLD*(1-exp(-MLD/h1));
        r = mumax*light/(40+light);
        
        N = N_max-y(1)-y(2);
        
        dydt = [
            %-(r*y(1)/(N0+y(1))-B(1))*y(2)+(1-a)*B(2)*y(3)*(y(2)^2/(y(2)^2+B(3)^2))+B(4)*y(3)^2;
            (r*N/(N0+N)-B(1))*y(1)-B(2)*y(2)*(y(1)^2/(y(1)^2+B(3)^2));
            a*B(2)*y(2)*(y(1)^2/(y(1)^2+B(3)^2))-B(4)*y(2)^2];
    end

end