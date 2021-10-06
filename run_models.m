function [pmlII,zmlII,pmlIII,zmlIII,dpdtII,dpdtIII,BII,BIII] = run_models
% function to solve biological models to produce figures in Freilich et al
% (in review) Biogeosciences. This functions solves for the concentration
% of zooplankton and phytoplankton using two different grazer function
% responses in a bulk mixed layer model

% Outputs:
% pmlII/pmlIII: mixed layer average phytoplankton concentration using a
% type II/type II function response
% zmlII/zmlIII: mixed layer average zooplankton concentration using a
% type II/type II function response
% dpdtII/dpdtIII: rate of change of phytoplankton biomass in the mixed
% layer using type II/type II function response
% BII/BIII: parameters used for the type II/type III model

%% parameters and options
% biological parameters in the order: [d_p g_0 p_0 d_z]
BIII = [0.001    4   15    1.8]; % parameters for type III model
BII = [0.0004    5.9   15    6]; % parameters for type II model
N_max = 30; % deep nutrient concentration
a = 0.5; % assimilation efficiency
mumax = 0.8; % maximum phytoplankton growth rate
N0 = 4; % nutrient half saturation parameter
x0 = [N_max 20 10]; % initial condition for ODE solver
options1 = odeset('Refine',1,'NonNegative',1); % options for ODE solver
t = linspace(0,3650,3650); % time interval over which to solve ODE
%% run models
[~,y1]=ode45(@f_typeII,t,x0,options1); % solve type II model
pmlII = y1(:,2)'; zmlII = y1(:,3)'; % save variables
[~,y1]=ode45(@f_typeIII,t,x0,options1); % sovle type III model
pmlIII = y1(:,2)'; zmlIII = y1(:,3)'; % save variables

% compute the change in biomass following Mignot et al. 2018 Nature
% Communications
[mld,~,~,~,~,~] = mldmodel(t);
dpdt1 = gradient(log(pmlII.*mld),mean(diff(t)));
dpdt2 = gradient(log(pmlII),mean(diff(t)));
dhdt = gradient(mld,mean(diff(t)));
dpdtII = dpdt1.*(dhdt > 0)+dpdt2.*(dhdt < 0);

[mld,~,~,~,~,~] = mldmodel(t);
dpdt1 = gradient(log(pmlIII.*mld),mean(diff(t)));
dpdt2 = gradient(log(pmlIII),mean(diff(t)));
dhdt = gradient(mld,mean(diff(t)));
dpdtIII = dpdt1.*(dhdt > 0)+dpdt2.*(dhdt < 0);

%% bulk mixed layer models

    function dydt = f_typeIII(t,y)
        h1 = 20;
        B = BIII;
        year_day=mod(t,365);
        loffset = 270;
        [MLD,~,x2,tm,ml_min,ml_max] = mldmodel(t);
        w_e=(ml_max-ml_min)*(-0.5*sin(x2)*pi/315/315*(tm-50)*2);
        w_e_pos=w_e;
        if w_e_pos<0
            w_e_pos=0;
        end
        
        yd1 = mod(year_day+loffset,365);
        h_light=20*(0.6*sin(yd1*pi/365*2)+1);
        light=h_light*h1/MLD*(1-exp(-MLD/h1));
        r = mumax*light/(40+light);
        
        dydt = [
            -(r*y(1)/(N0+y(1))-B(1))*y(2)+(1-a)*B(2)*y(3)*(y(2)^2/(y(2)^2+B(3)^2))+B(4)*y(3)^2+w_e_pos*(N_max-y(1))/MLD;
            (r*y(1)/(N0+y(1))-B(1))*y(2)-B(2)*y(3)*(y(2)^2/(y(2)^2+B(3)^2))-w_e_pos*y(2)/MLD;
            a*B(2)*y(3)*(y(2)^2/(y(2)^2+B(3)^2))-B(4)*y(3)^2-w_e_pos*y(3)/MLD];
    end

    function dydt = f_typeII(t,y)
        h1 = 20;
        B = BII;
        year_day=mod(t,365);
        loffset = 270;
        [MLD,~,x2,tm,ml_min,ml_max] = mldmodel(t);
        w_e=(ml_max-ml_min)*(-0.5*sin(x2)*pi/315/315*(tm-50)*2);
        w_e_pos=w_e;
        if w_e_pos<0
            w_e_pos=0;
        end
        
        yd1 = mod(year_day+loffset,365);
        h_light=20*(0.6*sin(yd1*pi/365*2)+1);
        light=h_light*h1/MLD*(1-exp(-MLD/h1));
        r = mumax*light/(40+light);
        
        dydt = [
            -(r*y(1)/(N0+y(1))-B(1))*y(2)+(1-a)*B(2)*y(3)*(y(2)/(y(2)+B(3)))+B(4)*y(3)^2+w_e_pos*(N_max-y(1))/MLD;
            (r*y(1)/(N0+y(1))-B(1))*y(2)-B(2)*y(3)*(y(2)/(y(2)+B(3)))-w_e_pos*y(2)/MLD;
            a*B(2)*y(3)*(y(2)/(y(2)+B(3)))-B(4)*y(3)^2-w_e_pos*y(3)/MLD];
    end

end