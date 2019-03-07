clear all; close all; clc;
%% Input parameters
% borehole
rb = [0.14/2, 0.15/2, 0.16/2, 0.17/2, 0.18/2. 0.19/2, 0.20/2]; % borehole radius [m]
db = rb * 2; % borehole diameter [m]
L = 1000; % borehole length [m]

% internal pipe
rpo = [0.07/2, 0.075/2, 0.08/2, 0.085/2, 0.09/2, 0.095/2, 0.10/2]; % internal pipe outer radius [m]
tp = [0.005, 0.006, 0.007, 0.008, 0.009, 0.01, 0.011]; % internal pipe thickness [m]
rpi = rpo-tp; % internal pipe inner raidus [m]
dpo = rpo*2; % internal pipe outer diameter [m]
kpp = 0.4; % internal pipe thermal conductivity [W/K-m] (PE)

% external pipe
reo = [0.115/2, 0.125/2, 0.135/2, 0.145/2, 0.155/2, 0.165/2, 0.175/2]; % external pipe outer radius [m]
te =  [0.005, 0.006, 0.007, 0.008, 0.009, 0.01, 0.011]; % external pipe thickness [m]
rei = reo-te; % external pipe inner raidus [m]
deo = 2*reo; % external pipe outer diameter [m]
kep = 54; % external pipe thermal conductivity [W/K-m] (Stainless)

% ground
kg = 2.6; % ground thermal conductivity [W/K-m]
cg = 2.35 * 10^6; % ground volumetric heat capacity [J/K-m3]
gg = 0.03; % geothermal gradient C/m

% grout/cement
kc = 2; % thermal conductivity of grout/cement [W/K-m]
cc = 3.35 * 10^6; % grout volumetric heat capacity [J/K-m3]

% fluid and flow
vfl = 2.5; % volumetric flow rate [l/s]
vf = vfl*0.001; % volumetric flow rate [m3/s]
rf = 977; % fluid density [kg/m3] (water)
cf = 4.19 * 10^6; % fluid volumetric heat capacity [J/K-m3]
kf = 0.65; % fluid thermal conductivity [W/K-m]
mf = 1.138 * 10^-3; % water viscosity [kg/m-s]
Pr = (cf/rf) * mf / kf; % Prandtl number for fluid

% Heat input rate
Q = [100000]; % heat input rate [W]

n_all=[];
COP_all=[];
dP_all=[];
Tbmean_all=[];

bh_info = [];

for i = 1:size(rb, 2)
    Tgs = 283.45; % ground srurface temp. [K] (10 degree celcius) - initial
    Qin = Q(1)
    rbin = rb(i)
    rpoin = rpo(i)
    reoin = reo(i)
    tpin = tp(i)
    tein= te(i)
    
    bh_info(i,:) = [rbin, rpoin, reoin, tpin, tein]
    
    Tgs = ((L* gg + Tgs)+Tgs)/2 % [K]
    [T1r, T2r, Tgr, Tcr, n , COPtotal, dP, Tbmean] = cxc_inj(rbin, L, rpoin, tpin, kpp, reoin, tein, kep, kg, Tgs, cg, gg, kc, cc, vfl,rf, cf, kf, mf, Qin);
    COP_all(i)= COPtotal
    n_all(i)= n
    dP_all(i)= dP
    Tbmean_all(i)= Tbmean
end

csvwrite(['COP.csv'], COP_all);
csvwrite(['efficiency.csv'], n_all);
csvwrite(['pressure_drop.csv'], dP_all);
csvwrite(['mean_Tb.csv'], Tbmean_all);
csvwrite(['bh_info.csv'], bh_info);


%% CXC----------------------------------
%% Assume constant heat injection rate 
function [T1r, T2r, Tgr, Tcr, n, COPtotal, dP, Tbmean] = cxc_inj(rb, L, rpo, tp, kpp, reo, te, kep, kg, Tgs, cg, gg, kc, cc, vfl,rf, cf, kf, mf, Q, tyr)

    %% Input parameters
    % borehole
    db = rb * 2; % borehole diameter [m]
    % internal pipe
    rpi = rpo-tp; % internal pipe inner raidus [m]
    dpo = rpo*2; % internal pipe outer diameter [m]
    % external pipe
    rei = reo-te; % external pipe inner raidus [m]
    deo = 2*reo; % external pipe outer diameter [m]
    % ground
    % grout/cement
    % fluid and flow
    vf = vfl*0.001; % volumetric flow rate [m3/s]
    Pr = (cf/rf) * mf / kf; % Prandtl number for fluid
    % Heat input rate
    %% Dimensionless Flow areas (1 = inflow, 2 = outflow)
    AD1 = rpi^2 / reo^2; % inflow area []
    AD2 = (rei^2-rpo^2)/(reo^2); % outflow area []

    %% Flow areas [m2]
    A1 = pi*rpi^2; % inflow, centre
    A2 = pi*rei^2 - pi*rpo^2;  % outflow, annulus

    %% Convective film coefficients
    hpi = cfc(2*rpi,Pr,kf,rf,mf,vfl,A1); % internal pipe convective film coefficient [W/K-m2]
    hei = cfc(2*rei-2*rpo,Pr,kf,rf,mf,vfl,A2); % external pipe convective film coefficient [W/K-m2]
    hpo = hei;

    %% Thermal resistances
    % thermal resistance between fluids in the internal pipe and the annulus
    R12 = 1/(2*pi*rpi*hpi) + 1/(2*pi*kpp) * log(rpo/rpi) + 1/(2*pi*rpo*hpo);
    % thermal resistance between fluids in the annulus and the grout
    Rg = 1/(2*pi*rei*hei) + 1/(2*pi*kep) * log(reo/rei);
    % thermal resistance of borehole 
    % (between the circulating fluid in the annulus and the borehole wall)
    Rb = 1/(2*pi*rei*hei) + 1/(2*pi*kep) * log(reo/rei) + 1/(2*pi*kc) * log(rb/reo);

    %% Other parameters
    Ns = (2*pi*kg*L)/(vf*cf); % dimensionless thermal conductance of ground
    Hf = cf/cg; % ratio of volumetric heat capacities of the circulating fluid
    N12 = L/(vf*cf*R12); % dimensionless thermal conductance of internal pipe
    Ng = L/(vf*cf*Rg); % dimensionless thermal conductance of grout
    Hg = cc/cg; % ratio of volumetric heat capacities of the grout and ground
    K = kc/kg; % ratio of the ground and ground thermal conductivities
    rDb = rb/reo;

    %% Analytical solution for CXC (injection in centre pipe)

    B1 = @(s) K.*sqrt(Hg.*s./K) .* besselk(1, rDb.*sqrt(Hg.*s./K)) .* besselk(0,rDb.*sqrt(s)) - sqrt(s).*besselk(0, rDb.*sqrt(Hg.*s./K)).*besselk(1,rDb.*sqrt(s));
    B2 = @(s) K.*sqrt(Hg.*s./K) .* besseli(1, rDb.*sqrt(Hg.*s./K)) .* besselk(0,rDb.*sqrt(s)) + sqrt(s).*besseli(0, rDb.*sqrt(Hg.*s./K)).*besselk(1,rDb.*sqrt(s));
    B4 = @(s) K.*sqrt(Hg.*s./K) .* besseli(1, rDb.*sqrt(Hg.*s./K)) .* besselk(0,rDb.*sqrt(Hg.*s./K)) + K.*sqrt(Hg.*s./K).*besseli(0, rDb.*sqrt(Hg.*s./K)).*besselk(1,rDb.*sqrt(Hg.*s/K));

    B2B1 = @(s) B2(s)/B1(s);
    B4B1 = @(s) B4(s)/B1(s);

    C0 = @(s) 1- 1./(1 - K.*Ns./Ng.* sqrt(Hg.*s./K).*(besseli(1, sqrt(Hg.*s./K)) - B2B1(s) .* besselk(1,sqrt(Hg.*s./K)))./(besseli(0, sqrt(Hg.*s./K)) + B2B1(s) .* besselk(0,sqrt(Hg.*s./K))));

    BB = @(s) (-Ns.*Hf.*AD2.*s)./2 - Ng.*C0(s) + (Ns.* Hf.* AD1.*s)./2;

    gam =@(s) -(Ns.* Hf.* AD2.*s./2 + N12 + Ng.*C0(s)) .* (Ns.*Hf.*AD1.*s./2+N12)+ N12.* N12;

    a1 = @(s) (-BB(s) + sqrt(BB(s).^2-4.*gam(s)))./2;
    a2 = @(s) (-BB(s) - sqrt(BB(s).^2-4.*gam(s)))./2;

    o1 = @(s) (a1(s) + Ns.*Hf.*AD1.*s./2 + N12)./N12;
    o2 = @(s) (a2(s) + Ns.*Hf.*AD1.*s./2 + N12)./N12;

    C1 = @(s) (Ns./s)./ ((1-o1(s)) .* (1-exp(a1(s)-a2(s))));
    C2 = @(s) ((Ns./s) - (1-o1(s)).*C1(s))./(1-o2(s));

    C3 = @(s) o1(s).*C1(s);
    C4 = @(s) o2(s).*C2(s);

    %% Time step and elevation step

    zs = [0:0.1:1]; % dimensionless depth ranging 0 to 1
    tts = [0:2620800:15724800]; % time step dt=1monnth
    ts = tts.*kg ./ (cg.*reo.^2); % converting the timestep to dimensioless 

    rrsg = [rb]; % radius from rb to 1 m for ground temperature distribution
    rsg = rrsg ./reo; % converting radius to dimensionless

    rrsc = [reo]; % radius from reo to rb for grout temperature distribution
    rsc = rrsc ./reo; % converting radius to dimensionless

    %% for loops
    T1 = zeros(size(zs,2),size(ts,2));
    T2 = zeros(size(zs,2),size(ts,2));
    Tca = [];
    Tga = [];

    for i = 1:size(zs,2)
        fprintf('looping...\n')
        zD = zs(i);
        TD1L = @(s) C1(s).*exp(a1(s).*zD) + C2(s).*exp(a2(s).*zD);
        TD2L = @(s) C3(s).*exp(a1(s).*zD) + C4(s).*exp(a2(s).*zD);
        B1zD = @(s) (Ng/(K*Ns) * TD1L(s)) / (Ng/(K*Ns)* (besseli(0,sqrt(Hg*s/K)) + B2B1(s)* besselk(0,sqrt(Hg*s/K)))-sqrt(Hg*s/K)*(besseli(1,sqrt(Hg*s/K)) - B2B1(s)* besselk(1,sqrt(Hg*s/K)))); 

        T1t = talbot_inversion(TD1L,ts)';
        T2t = talbot_inversion(TD2L,ts)'; 
        T1(i,:) = T1t; % T1 = annulus dimensionless temperature
        T2(i,:) = T2t; % T2 = centre dimensionless temperature

        for j=1:size(rsg,2)
            rDg = rsg(j);
            TgL = @(s) B1zD(s) * B4B1(s) * besselk(0,rDg*sqrt(s));
            Tgt = talbot_inversion(TgL,ts)';
            Tga(i,:,j) = Tgt; % Tga = ground dimensionless temperature

            for n = 1:size(rsc,2)
                rDc = rsc(n);
                TcL = @(s) B1zD(s) * (besseli(0,rDc*sqrt(Hg*s/K)) + B2B1(s)* besselk(0,rDc*sqrt(Hg*s/K)));
                Tct = talbot_inversion(TcL,ts)'; 
                Tca(i,:,n) = Tct; % Tca = grout dimensionless temperature
            end 
        end
    end

    %% Converting dimensionless unit to appropriate unit for plotting
    rgr = rrsg; % [m]
    rcr = rrsc; % [m]
    tr = tts .* 1.1574e-5; % [day]
    zr = zs.*L; % [m]
    Tg = Tgs; % [K](Undisturbed mean ground temperature)

    for i = 1:size(T1,2)
        T1r(:,i) = (T1(:,i).*Q)./(2.*pi.*kg.*L) + Tg - 273.15; %[C]
        T2r(:,i) = (T2(:,i).*Q)./(2.*pi.*kg.*L) + Tg - 273.15; %[C]
        for j = 1:size(Tga,3)
            Tgr(:,i,j) = (Tga(:,i,j).*Q)./(2.*pi.*kg.*L) + Tg - 273.15; %[C]
            Tcr(:,i,j) = (Tca(:,i,j).*Q)./(2.*pi.*kg.*L) + Tg - 273.15; %[C]
        end
    end
    
    n = 1 - 10/mean(Tgr(:,end))
    [COPtotal, dP] = COPtot(2*rei-2*rpo,rf,mf,vfl,A2,L,Q);
    COPtotal = COPtotal
    dP= dP
    Tbmean = mean(Tgr(:,end))

    %% save matrices to files
    csvwrite(['time_cxc_' num2str(rb) '.csv'] ,tr);
    csvwrite(['centre_pipe_in_temp_cxc_' num2str(rb) '.csv'], T1r);
    csvwrite(['annulus_pipe_out_temp_cxc_' num2str(rb) '.csv'], T2r);
    csvwrite(['ground_temp_cxc_' num2str(rb) '.csv'], Tgr);
    csvwrite(['grout_temp_cxc_' num2str(rb) '.csv'], Tcr);
    csvwrite(['depth_cxc_' num2str(rb) '.csv'], zr);
    csvwrite(['radius_ground_cxc_' num2str(rb) '.csv'], rgr);
    csvwrite(['radius_grout_cxc_' num2str(rb) '.csv'], rcr);

    %% Functions
    function h = cfc(d,Pr,kf,rf,mf,Vf,A)
    % Function determining convective film coefficients, h [W/K-m2]
        V = 0.001 * Vf / A; % converting l/s to m3/s to m/s    
        Re = rf*V*d/mf;

        if Re < 2300
            Nu = 4.364;
            h = Nu * kf /d;
        else
            % f = friction factor
            f = 2 / ((((8/Re)^10 + (Re/36500)^20)^-0.5 + (2.21*log(Re/7))^10)^(1/5));
            Nu = ((f/2)*(Re-1000)*Pr)/(1+12.7*(f/2)^0.5 * (Pr^(2/3) - 1));
            h = Nu * kf /d;
        end
    end

    function [COP, dP] = COPtot(d,rf,mf,Vf,A,l,Q)
    % Function determining convective film coefficients, h [W/K-m2]
        V = 0.001 * Vf / A; % converting l/s to m3/s to m/s
        mfl = 0.001*Vf*rf;
        Dh = d;
        Re = rf*V*d/mf;
        
        f = (0.790*log(Re)-1.64)^(-2);
        dP = f*l /Dh * rf *V^2/2;
        Wpump = dP * mfl / (rf*0.75);
        Whp = Q/(4-1);
        COP = (Q + Whp) / (Whp + Wpump);
        dP = dP*(1e-5);
    end

    function ilt = talbot_inversion(f_s, t, M)
    % ilt = talbot_inversion(f_s, t, [M])
    %
    % Returns an approximation to the inverse Laplace transform of function
    % handle f_s evaluated at each value in t (1xn) using Talbot's method as
    % summarized in the source below.
    % 
    % This implementation is very coarse; use talbot_inversion_sym for better
    % precision. Further, please see example_inversions.m for discussion.
    %
    % f_s: Handle to function of s
    % t:   Times at which to evaluate the inverse Laplace transformation of f_s
    % M:   Optional, number of terms to sum for each t (64 is a good guess);
    %      highly oscillatory functions require higher M, but this can grow
    %      unstable; see test_talbot.m for an example of stability.
    % 
    % Abate, Joseph, and Ward Whitt. "A Unified Framework for Numerically 
    % Inverting Laplace Transforms." INFORMS Journal of Computing, vol. 18.4 
    % (2006): 408-421. Print.
    % 
    % The paper is also online: http://www.columbia.edu/~ww2040/allpapers.html.
    % 
    % Tucker McClure
    % Copyright 2012, The MathWorks, Inc.
        % Make sure t is n-by-1.
        if size(t, 1) == 1
            t = t';
        elseif size(t, 2) > 1
            error('Input times, t, must be a vector.');
        end
        % Set M to 64 if user didn't specify an M.
        if nargin < 3
            M = 64;
        end

        % Vectorized Talbot's algorithm

        k = 1:(M-1); % Iteration index

        % Calculate delta for every index.
        delta = zeros(1, M);
        delta(1) = 2*M/5;
        delta(2:end) = 2*pi/5 * k .* (cot(pi/M*k)+1i);

        % Calculate gamma for every index.
        gamma = zeros(1, M);
        gamma(1) = 0.5*exp(delta(1));
        gamma(2:end) =    (1 + 1i*pi/M*k.*(1+cot(pi/M*k).^2)-1i*cot(pi/M*k))...
                       .* exp(delta(2:end));

        % Make a mesh so we can do this entire calculation across all k for all
        % given times without a single loop (it's faster this way).
        [delta_mesh, t_mesh] = meshgrid(delta, t);
        gamma_mesh = meshgrid(gamma, t);

        % Finally, calculate the inverse Laplace transform for each given time.
        ilt = 0.4./t .* sum(real(   gamma_mesh ...
                                 .* arrayfun(f_s, delta_mesh./t_mesh)), 2);
    end

    function ilt = euler_inversion(f_s, t, M)
    % ilt = euler_inversion(f_s, t, [M])
    %
    % Returns an approximation to the inverse Laplace transform of function
    % handle f_s evaluated at each value in t (1xn) using the Euler method as
    % summarized in the source below.
    % 
    % This implementation is very coarse; use euler_inversion_sym for better
    % precision. Further, please see example_inversions.m for examples.
    %
    % f_s: Handle to function of s
    % t:   Times at which to evaluate the inverse Laplace transformation of f_s
    % M:   Optional, number of terms to sum for each t (64 is a good guess);
    %      highly oscillatory functions require higher M, but this can grow
    %      unstable; see test_talbot.m for an example of stability.
    % 
    % Abate, Joseph, and Ward Whitt. "A Unified Framework for Numerically 
    % Inverting Laplace Transforms." INFORMS Journal of Computing, vol. 18.4 
    % (2006): 408-421. Print.
    % 
    % The paper is also online: http://www.columbia.edu/~ww2040/allpapers.html.
    % 
    % Tucker McClure
    % Copyright 2012, The MathWorks, Inc.
        % Make sure t is n-by-1.
        if size(t, 1) == 1
            t = t';
        elseif size(t, 2) > 1
            error('Input times, t, must be a vector.');
        end
        % Set M to 64 if user didn't specify an M.
        if nargin < 3
            M = 32;
        end

        % Vectorized Talbot's algorithm
        bnml = @(n, z) prod((n-(z-(1:z)))./(1:z));

        xi = [0.5, ones(1, M), zeros(1, M-1), 2^-M];
        for k = 1:M-1
            xi(2*M-k + 1) = xi(2*M-k + 2) + 2^-M * bnml(M, k);
        end
        k = 0:2*M; % Iteration index
        beta = M*log(10)/3 + 1i*pi*k;
        eta  = (1-mod(k, 2)*2) .* xi;

        % Make a mesh so we can do this entire calculation across all k for all
        % given times without a single loop (it's faster this way).
        [beta_mesh, t_mesh] = meshgrid(beta, t);
        eta_mesh = meshgrid(eta, t);

        % Finally, calculate the inverse Laplace transform for each given time.
        ilt =    10^(M/3)./t ...
              .* sum(eta_mesh .* real(arrayfun(f_s, beta_mesh./t_mesh)), 2);
    end

end


