Nt = 1200; % Number of time points
Ns = 4; % Number of state variables
No = Ns-1; % Number of observable state variables
Np = 6; % Number of parameters

%------ TIME GRID ------%
tmin = 0;
tmax = 360;
t = linspace(tmin, tmax, Nt);
h = (tmax - tmin)/Nt;

%------ PARAMETER ESTIMATES ------%
% psid psir gamma ed er KL
p=[0.0392, 0.0571, 1.27e-10, 1.86e-9, 6.22e-10, 145];

%------ INITIAL CONDITIONS ------%
% recipients donors transconjugants substrates
r0 = 1.96*1e6;
d0 = 1.77*1e6;
t0 = 0.0;
s0 = 1.0;

% state vector
x = zeros(Ns,Nt); 
x(:,1) = [r0, d0, t0, s0];

%------ CALCULATE STATE VECTOR ------%
x = experiment(x,p);

%------ CALCULATE LOCAL SENSITIVITIES ------%
xs = zeros(Ns,Nt,Np);
for pn=1:Np
    fp = ones(1,Np);
    fp(pn) = fp(pn) + 0.01; % forward step vector
    xf = zeros(Ns,Nt); xf(:,1) = [r0, d0, t0, s0];% forward step state vector
    xf = experiment(xf,fp.*p); % forward step in parameter pn
    xs(:,:,pn) = (xf./x-1)/0.01; % finite difference derivative for normalized sensitivities
end
xs = permute(xs,[1,3,2]);

%------ OUTPUT SENSITIVITIES ------%
% format sensitivities for output
xso = zeros(Ns*Np, Nt);
for k=1:Nt
    for p=1:Np
        for s=1:Ns
            if xs(s,p,k) < 1e8 % doesn't include x=0
                xso(s + Ns*(p-1), k) = xs(s,p,k);
            end
        end
    end
end
% output sensitivities
dlmwrite('Sensitivities/sensitivities1.txt', xso, 'delimiter', ' ', 'precision', 6)

%------ PLOTTING ------%
plottingStates(x);
%plottingSensitivities(xs);
%plottingSquaredSensitivities(xs);
%plottingFIMEigenvalues(xs);

%------ INTEGRATED SQUARED SENSITIVITIES ------%
S = zeros(1,Np);
for k = 1:Nt
    for p = 1:Np
        for s = 1:No
            if xs(s,p,k) < 1e8 % doesn't include sensitivites where x=0
                S(p) = S(p) + xs(s,p,k)^2;
            end
        end
    end
end
disp("Sensitivites are")
disp(S/Nt);

%====== FUNCTIONS ======%

%------ EXPERIMENT ------%
function x = experiment(x,p)
    %------ TIME GRID ------%
    tmin = 0;
    tmax = 360;
    Nt = 1200;
    h = (tmax - tmin)/Nt;

    %------ PARAMETERS ------%
    psid = p(1);    psir = p(2);    % max growth rates
    gamma = p(3);                   % max cg rate
    ed = p(4);      er = p(5);      % yield
    KL = p(6);                      % growth lag time

    %------ EULER'S METHOD ------%
    for tn = 1:(Nt-1)
        ti = tmin + h*tn;
        % state variables
        rp = x(1,tn);  
        dn = x(2,tn);
        tc = x(3,tn);  
        sb = x(4,tn);
    
        % growth factors
        dn_g = psid*(ti/(KL+ti))*(sb);
        rp_g = psir*(ti/(KL+ti))*(sb);
    
        tc_g = 0;
        cg_rate = gamma*sb*(dn+tc)*rp;
    
        % model rhs
        dxdt = [rp_g*rp - cg_rate;
                dn_g*dn;
                tc_g*tc + cg_rate;               
                -er*rp_g*rp - ed*dn_g*dn - er*tc_g*tc];
        
        x(:,tn+1) = x(:,tn) + h*dxdt; 
    end
end

%------ PLOTTING ------%
%#ok<*DEFNU>

function plottingStates(x)

    Nt = 1200; % Number of time points

    %------ PLOTTING ------%
    figure(1)
    clf;
    hold on

    %------ TIME GRID ------%
    tmin = 0;
    tmax = 360;
    t = linspace(tmin, tmax, Nt);
    
    % plotting state variables
    plot(t, log10(x(1,:)), '-', 'Color', 'k','LineWidth',2)
    plot(t, log10(x(2,:)), '-', 'Color', [36, 7, 133]/256,'LineWidth',2)
    plot(t, log10(x(3,:)), '-', 'Color', [168, 26, 0]/256,'LineWidth',2)
    ax = gca;
    ax.FontSize = 24; 
    %xlabel('Time (min)')
    xlabel(' ')
    ylabel('Abundance (cells/cm^2)')
    title('Optimal Experiment #1','FontSize',24)
    axis([tmin-12 tmax 4 9]);
    xticks([0 100 200 300])
    xticklabels({'0' '100' '200' '300'});
    yticks([4 6 8])
    yticklabels({'10^4' '10^6' '10^8'})
    %legend({'recipients','donors','transconjugants'},'Location','southeast','FontSize',24)

    % plotting optimal sampling intervals (shaded region)
    si_opt = [16, 26, 34, 35, 39, 40]; % Optimal sampling intervals (shaded region)
    for mi=1:6
        x1 = linspace(9*(si_opt(mi)-1),9*si_opt(mi),20);
        x2 = 10*(x1+1)./(x1+1);
        fill([x1 flip(x1)],[x2 zeros(size(x2))],[36, 7, 133]/256,'LineStyle','none')
        plot(x1,x2,'k-');
        alpha(0.25);
    end
    legend({'recipients','donors','transconjugants'},'Location','southeast','FontSize',24)
    box on
end
    
function plottingSensitivities(xs)
    Nt = 1200; % Number of time points
    Np = 6; % Number of parameters
    
    %------ TIME GRID ------%
    tmin = 0;
    tmax = 360;
    t = linspace(tmin, tmax, Nt);
    
    %------ PLOTTING SENSITIVITIES ------%
    p_str = ["\psi^{max}_D", "\psi^{max}_R", "\gamma^{max}", "\epsilon_D", "\epsilon_R", "K_L"];
    for p=1:Np
        figure(p+1); clf; hold on;
        plot(t, xs(1,:,p), '-', 'Color', 'k','LineWidth',2)
        plot(t, xs(2,:,p), '-', 'Color', [36, 7, 133]/256,'LineWidth',2)
        plot(t, xs(3,:,p), '-', 'Color', [168, 26, 0]/256,'LineWidth',2)
    
        si = [16, 26, 34, 35, 39, 40];
        for mi=1:6
            x1 = linspace(9*(si(mi)-1),9*si(mi),20);
            x2 = 20*(x1+1)./(x1+1);
            axis([0 360 -4 6]);
            fill([x1 flip(x1)],[x2 -20*ones(size(x2))],[36, 7, 133]/256,'LineStyle','none')
            plot(x1,x2,'k-');
            alpha(0.25);
        end
    
        ax = gca;
        ax.FontSize = 32; 
        %xlabel('Time (min)')
        xlabel('')
        %ylabel('Sensitivity')
        ylabel('')
        xticks([0 100 200 300])
        xticklabels({'0' '100' '200' '300'});
        yticks([-4 -2 0 2 4 6])
        yticklabels({'-4' '-2' '0' '2' '4' '6'})
        %ylims = ylim;
        text(10,4.5, p_str(p), 'FontSize',64);
        box on
        
        figure(2);
        legend({'recipients','donors','transconjugants'},'Location','southwest','FontSize',32)
    end
end

function plottingSquaredSensitivities(xs) 
    Nt = 1200; % Number of time points
    Np = 6; % Number of parameters
    
    %------ TIME GRID ------%
    tmin = 0;
    tmax = 360;
    t = linspace(tmin, tmax, Nt);
    
    %------ PLOTTING SQUARED SENSITIVITIES ------%
    p_str = ["\psi^{max}_D", "\psi^{max}_R", "\gamma^{max}", "\epsilon_D", "\epsilon_R", "K_L"];
    for p=1:Np
        figure(p+1); clf; hold on;
        plot(t, xs(1,:,p).*xs(1,:,p), '-', 'Color', 'k','LineWidth',2)
        plot(t, xs(2,:,p).*xs(2,:,p), '-', 'Color', [36, 7, 133]/256,'LineWidth',2)
        plot(t, xs(3,:,p).*xs(3,:,p), '-', 'Color', [168, 26, 0]/256,'LineWidth',2)
    
        si = [16, 26, 34, 35, 39, 40];
        for mi=1:6
            x1 = linspace(9*(si(mi)-1),9*si(mi),30);
            x2 = 40*(x1+1)./(x1+1);
            fill([x1 flip(x1)],[x2 -20*ones(size(x2))],[36, 7, 133]/256,'LineStyle','none')
            plot(x1,x2,'k-');
            alpha(0.25);
        end
    
        axis([0 360 0 30]);
        ax = gca;
        ax.FontSize = 32; 
        %xlabel('Time (min)')
        xlabel('')
        %ylabel('Sensitivity')
        ylabel('')
        xticks([0 100 200 300])
        xticklabels({'0' '100' '200' '300'});
        yticks([0 6 12 18 24 30])
        yticklabels({'0' '6' '12' '18' '24' '30'})
        %legend({'recipients','donors','transconjugants'},'Location','southwest','FontSize',32)
        %ylims = ylim;
        text(10,6*4.5, p_str(p), 'FontSize',64);
        box on
    end
end

function plottingFIMEigenvalues(xs) 
    Nt = 1200; % Number of time points
    Np = 6; % Number of parameters
    No = 3; % Number of observable state variables
    
    %------ TIME GRID ------%
    tmin = 0;
    tmax = 360;
    t = linspace(tmin, tmax, Nt);
    h = (tmax - tmin)/Nt;
    
    %------ CALCULATE FIM EIGENVALUES ------%
    sv = zeros(Np,No);
    lam = zeros(No,Nt);
    FIM = zeros(Np,Np);
    for ti=2:Nt
        for p=1:Np
            for s=1:No
                if xs(s,p,ti) < 1e8
                    sv(s,p) = xs(s,p,ti);
                end
            end
        end
 
        FIM = FIM + h*(sv'*sv);
        [~,D] = eig(sv'*sv);
        for p=1:Np 
            lam(p,ti) = D(p,p); 
        end
    end

    
    %[V,D] = eig(sv(:,700)*sv(:,700)');
    
    %------ PLOTTING FIM EIGENVALUES ------%
    p_str = ["\psi^{max}_D", "\psi^{max}_R", "\gamma^{max}", "\epsilon_D", "\epsilon_R", "K_L"];
    for p=1:Np
        figure(p+1); clf; hold on;
        plot(t, log(lam(p,:)), '-', 'Color', 'k','LineWidth',2)
        %plot(t, log(lam(2,:)), '-', 'Color', [36, 7, 133]/256,'LineWidth',2)
        %plot(t, log(lam(3,:)), '-', 'Color', [168, 26, 0]/256,'LineWidth',2)
    
        si = [16, 26, 34, 35, 39, 40];
        for mi=1:6
            x1 = linspace(9*(si(mi)-1),9*si(mi),30);
            x2 = 40*(x1+1)./(x1+1);
            fill([x1 flip(x1)],[x2 -20*ones(size(x2))],[36, 7, 133]/256,'LineStyle','none')
            plot(x1,x2,'k-');
            alpha(0.25);
        end
    
        %axis([0 360 0 30]);
        ax = gca;
        ax.FontSize = 32; 
        %xlabel('Time (min)')
        xlabel('')
        %ylabel('Sensitivity')
        ylabel('')
        %xticks([0 100 200 300])
        %xticklabels({'0' '100' '200' '300'});
        %yticks([0 6 12 18 24 30])
        %yticklabels({'0' '6' '12' '18' '24' '30'})
        %legend({'recipients','donors','transconjugants'},'Location','southwest','FontSize',32)
        %ylims = ylim;
        text(10,24, p_str(p), 'FontSize',64);
        box on
    end
    
    figure(2)
    %legend({'recipients','donors','transconjugants'},'Location','southwest','FontSize',32)
end
