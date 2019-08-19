tic

Nt = 1200;  % Number of time points
Ns = 4;     % Number of state variables
No = Ns-1;  % Number of observable state variables
Np = 6;     % Number of parameters
Nm = 6;     % Number of measurements
Nk = 40;    % Number of measurement intervals (for now: must be a factor of Nt)

%------ LOAD SENSITIVITIES ------%
S = load('Sensitivities/sensitivities1.txt');
S = reshape(S,Ns,Np,Nt);
Sk = zeros(No, Np, Nk); % average sensitivities for all Nk intervals
for k=1:Nk
    iL = ceil(Nt/Nk); % interval length
    Sk(:,:,k) = sum( S(1:No, :, iL*(k-1)+1:iL*k), 3)/iL;
end

%------ CALCULATING J ------%
C = combnk(1:Nk,Nm);
Nc = size(C,1); % Number of configurations
J = zeros(1,Nc);
for c=1:Nc
    mS = Sk(:,:,C(c,:));
    FIM = zeros(Np,Np);
    for mi=1:Nm
       FIM = FIM + mS(:,:,mi)'*mS(:,:,mi);
    end
    J(c) = det(FIM);
end

%------ SORTING ------%
[Js,I] = sort(J, 'descend');

toc

disp(strcat('c*=', num2str(I(1))));

%------ OUTPUT ------%
outfile = "Optimality/" + "J_Nk=" + num2str(Nk) + ".txt";
dlmwrite(outfile, [Js;1:Nc] , 'delimiter', ' ', 'precision', 6);

%------ PLOTTING ------%
%plottingAll(Js,I,Nc);

%====== FUNCTIONS ======%
%{
%------ PLOTTING ------%
function plottingAll(Js,I,Nc)

    figure(1)

    clf;
    hold on
    %cgrid = 1:Nc;
    %cgrid2 = 10.^(linspace(0,7,50));
    %A = 445;
    %b = 0.35;
    %c = 60000;
    %strexp_fit = A*exp(-(log(cgrid)/c));
    %strexp_fit = A*exp(-((cgrid2-1)/c).^b);
    %{
    A20 = J20(1); cr20 = 0.037*size(J20,2); b20 = 0.48;
    fcgrid20 = 10.^(linspace(0,log10(size(J20,2)),50));
    strexp_fit20 = A20*exp(-((fcgrid20-1)/cr20).^b20);
    plot(fcgrid20, strexp_fit20, 'o', 'Color', [1, 0, 98]/256,'LineWidth',1.5);
    plot(cgrid20, J20, '-', 'Color', [1, 0, 98]/256,'LineWidth',2);
    
    A30 = J30(1); cr30 = 0.023*size(J30,2); b30 = 0.45;
    fcgrid30 = 10.^(linspace(0,log10(size(J30,2)),50));
    strexp_fit30 = A30*exp(-((fcgrid30-1)/cr30).^b30);
    plot(fcgrid30, strexp_fit30, 'o', 'Color', [0, 106, 110]/256,'LineWidth',1.5);
    plot(cgrid30, J30, '-', 'Color', [0, 106, 110]/256,'LineWidth',2);
    
    A40 = J40(1); cr40 = 56000; b40 = 0.35;
    fcgrid40 = 10.^(linspace(0,log10(size(J40,2)),50));
    strexp_fit40 = A40*exp(-((fcgrid40-1)/cr40).^b40);
    plot(fcgrid40, strexp_fit40, 'o', 'Color', [2, 120, 0]/256,'LineWidth',1.5);
    plot(cgrid40, J40, '-', 'Color', [2, 120, 0]/256,'LineWidth',2);
    
    A50 = J50(1); cr50 = 180000; b50 = 0.35;
    fcgrid50 = 10.^(linspace(0,log10(size(J50,2)),50));
    strexp_fit50 = A50*exp(-((fcgrid50-1)/cr50).^b50);
    plot(fcgrid50, strexp_fit50, 'o', 'Color', [127, 74, 0]/256,'LineWidth',1.5);
    %plot(cgrid50, J50, '-', 'Color', [127, 74, 0]/256,'LineWidth',2);
%}
    A60 = J60(1); cr60 = 472000; b60 = 0.35;
    fcgrid60 = 10.^(linspace(0,log10(size(J60,2)),50));
    strexp_fit60 = A60*exp(-((fcgrid60-1)/cr60).^b60);
    plot(fcgrid60, strexp_fit60, 'o', 'Color', [130, 0, 0]/256,'LineWidth',1.5);
    plot(cgrid60, J60, '-', 'Color', [130, 0, 0]/256,'LineWidth',2);

    %plot(cgrid, Js, '-', 'Color', [130, 0, 0]/256,'LineWidth',2);
    %plot(cgrid, Js, '-', 'Color', [131, 0, 73]/256,'LineWidth',2);
    %plot(cgrid2, strexp_fit, 'o', 'Color', 'k','LineWidth',1.5);
    ax = gca;
    ax.FontSize = 24; 
    xlabel('sampling schedule configuration','FontSize',24)
    ylabel('J(c)','FontSize',24)
    set(gca, 'XScale', 'log');
    axis([1 1e8 0 600]);
    xticks([1e0 1e2 1e4 1e6 1e8 1e10])
    xticklabels({'10^0' '10^2' '10^4' '10^6' '10^8' '10^10'})
    yticks([0 200 400 600])
    yticklabels({'0' '200' '400' '600'})
    %ylims = ylim;
    %text(2,60, strcat('c*=', num2str(I(1))), 'FontSize',32, 'fontweight','bold');
    %text(2,30, strcat('J(c*)=', num2str(253.4)), 'FontSize',32, 'fontweight','bold');
    %legend('N_k=20', 'N_k=30', 'N_k=40', 'N_k=50', 'N_k=60','Location','northeast')
    box on
end
%}
% blue [36, 7, 133]/256