%% ########### Evaluation of estimation results for UKF-(AW)GSI and D-UKF-AW-GSI at L-TOWN network ########## %%

clear; clc; rng('default')
addpath('functions');

%% Configuration %%

networkID = 'Modena'; %'LTOWN';
hyd_data_file = 'leaktionary_Modena_45to65_1p100_testing.mat';%'estimation_required_data.mat';%'2018_SCADA_Measurements_zoneA_withdemands.mat';
struc_data_file = {'ModenaNewGraphInfo.mat',...
                   'ModenaNewIncidenceMatrix.mat',...
                   'ModenaNewPipeDistance.mat',...
                   'sensors.mat'};%{'coordenades_BattleDIM1_1.mat'};
enable_plots = [1 1 1]; cnt_plots = 1;

% Load datasets %

load(fullfile('data',networkID,'hydraulic',hyd_data_file));
for i=1:length(struc_data_file)
    load(fullfile('data',networkID,'structure',struc_data_file{i}));
end

% Structural data %

N = length(A);

%%% (Only for L-TOWN) %%%

if strcmp(networkID,'LTOWN')

    reservoirsID = [111 300];%[d.NodeReservoirIndex d.NodeTankIndex];

    zone_str = 'zoneA';
    reservoirsID_comp = reservoirsID;
    N_comp = N;

    % Divide complete node set in the different zones %

    zoneA = []; zoneB = []; zoneC = [];

    for i=1:N_comp
        if d.NodeElevations(i)<16
            zoneB = [zoneB;i];
        elseif d.NodeElevations(i)>=16 && d.NodeElevations(i)<32
            zoneA = [zoneA;i];
        elseif d.NodeElevations(i)>=32 && d.NodeElevations(i)<48
            zoneA = [zoneA;i];
        elseif d.NodeElevations(i)>=48 && d.NodeElevations(i)<64
            zoneC = [zoneC;i];
        else
            zoneC = [zoneC;i];
        end
    end

    zoneAorig = zoneA;
    zoneBorig = zoneB;
    zoneCorig = zoneC;

    zoneA = [zoneA;783;784];

    toDeleteNodesFromGraph = [d.NodeReservoirIndex 336 303];

    for i=1:length(toDeleteNodesFromGraph)
        fA = find(zoneA == toDeleteNodesFromGraph(i));
        zoneA(fA) = [];
        fB = find(zoneB == toDeleteNodesFromGraph(i));
        zoneB(fB) = [];
        fC = find(zoneC == toDeleteNodesFromGraph(i));
        zoneC(fC) = [];
    end

    if strcmp(zone_str,'zoneA')
        zone_to_study = zoneA;
        disp('Zone A')
    elseif strcmp(zone_str,'zoneB')
        zone_to_study = zoneB;
        disp('Zone B')
    else
        zone_to_study = zoneC;
        disp('Zone C')
    end

    % Extract pipes of Zone %

    E_comp = d.NodesConnectingLinksIndex;

    pipes_in_zoneA = []; Ec = []; k=1;
    for pipe=1:length(E_comp)
        if ~isempty(intersect(zone_to_study,E_comp(pipe,1))) && ~isempty(intersect(zone_to_study,E_comp(pipe,2)))
            pipes_in_zoneA = [pipes_in_zoneA;pipe];
            E(k,:) = [find(zone_to_study==E_comp(pipe,1)) find(zone_to_study==E_comp(pipe,2))];
            k=k+1;
        end
    end

    nEdges = length(E);

    % Derive parameters associated to Zone A %

    A = A(zone_to_study,zone_to_study);
    N = length(zone_to_study);
    node_coordenades = node_coordenades(zone_to_study,:);

    reservoirsID_A = intersect(zone_to_study,reservoirsID_comp);
    for i=1:length(reservoirsID_A)
        reservoirsID(1,i) = find(zone_to_study==reservoirsID_A(i));
    end
    PipeDistance = PipeDistance(zone_to_study,zone_to_study);

    % Sensors %

    sensors_comp = [sensors;reservoirsID_A]; % add the tank
    sensors_in_zoneA = intersect(zoneA,sensors_comp);
    sensors_a = zeros(length(sensors_in_zoneA),1);
    for s=1:length(sensors_in_zoneA)
        sensors_a(s) = find(zoneA==sensors_in_zoneA(s));
    end
    sensors = sensors_a;

    % Create graph %

    G = graph(A);
    GEdges = table2array(G.Edges); GEdges = GEdges(:,1:2); E_orig=E;
    for edge=1:length(GEdges)
        G2E = find(ismember(E,[GEdges(edge,1) GEdges(edge,2)],'rows') + ismember(E,[GEdges(edge,2) GEdges(edge,1)],'rows')==1);
        G.Edges.Weight(edge) = d.LinkLength(G2E);
    end

    % Derive the weighted adjacency matrix %

    WA = (PipeDistance.*A);
    WA=1./WA;
    WA(WA==Inf)=0;

    % Derive incidence matrix (aproximated by topology only) %

    innernodesIDs = setdiff(1:N,reservoirsID);
    C = zeros(size(A));
    for i=1:length(reservoirsID)
        for j = 1:length(innernodesIDs)
            path = shortestpath(G,reservoirsID(i),innernodesIDs(j));
            for k = 1:length(path)-1
                C(path(k),path(k+1)) = C(path(k),path(k+1)) + 1;
            end
        end
    end
    I = zeros(nEdges,N);
    for i=1:length(E)
        if C(E(i,1),E(i,2)) > C(E(i,2),E(i,1))
            I(i,E(i,1)) = -1; I(i,E(i,2)) = 1;
        else
            I(i,E(i,2)) = -1; I(i,E(i,1)) = 1;
        end
    end

elseif strcmp(networkID,'Modena')
    node_coordenades = nc;
    time_instant = 1; % this selects a leak of 4.5 l/s
    % Prepare data %
    analysed_leaks = 1:N-length(reservoirsID);
    for j=1:length(analysed_leaks)
        hnom_M(:,j) = Leaktionary_nom{j}.Head(time_instant,:)';
        hleak_M(:,j) = Leaktionary{j}.Head(time_instant,:)';
        qleak_M(:,j) = Leaktionary{j}.Flows(time_instant,:)'/1000; % l/s a m3/s
        dleak_M(:,j) = Leaktionary{j}.Demand(time_instant,:)'/1000; % l/s a m3/s
    end
else
    warning('If you are setting a custom network, check the data structures.')
end

% (optional) Plot network %

if enable_plots(cnt_plots)==1
    figure(cnt_plots);
    fig = plot(G,'XData',node_coordenades(:,1),'YData',node_coordenades(:,2),'Marker','o','MarkerSize',4);
    xlabel('Latitude','interpreter','latex','fontsize',14);
    ylabel('Longitude','interpreter','latex','fontsize',14);
    % title(sprintf('%s - Network graph',network_name));
    grid

    highlight(fig,reservoirsID,'MarkerSize',7,'NodeColor','g','Marker','s');
    highlight(fig,sensors,'MarkerSize',7,'NodeColor','m');%,'Marker','s');
end
cnt_plots = cnt_plots + 1;

%% Graph-based State Interpolation - Configuration %%

% Settings %

tau = 1000;

% Generate sensor location matrix %

s = zeros(N,1); s(sensors) = 1;
S = eye(N).*(s*s');
Ss = [S(sensors,:)];

% Generate Laplacian matrix %

WA = sparse(WA);
Deg = diag(sum(WA));
Lap2 = Deg - WA;
Ld = Lap2'*Deg^-2*Lap2;

% Generate resistance coefficient matrix %

lengths = d.LinkLength; 
diameters = d.LinkDiameter/1000; 
roughness = d.LinkRoughnessCoeff; 

if strcmp(networkID,'LTOWN') % to correct some issues with the modelling affecting GSI
    lengths(end-3:end) = 1; roughness(end-3:end) = 140;
    pipes_list = pipes_in_zoneA;
elseif strcmp(networkID,'Modena')
    pipes_list = 1:length(E);
end 

t = 10.674*lengths./(diameters.^4.871.*roughness.^1.852);
T = diag(t(pipes_list)); invT = inv(T);

%% Estimation %%

if strcmp(networkID,'LTOWN')
    vsensors_file = 'virtual_sensors_100(657nodes).mat';
elseif strcmp(networkID,'Modena')
    vsensors_file = 'virtual_sensors_40.mat'; 
end
load(fullfile('data',networkID,'structure',vsensors_file));
idk = vsensors;

for index = 1:length(analysed_leaks)
    disp(['######## Leak ' num2str(index) ' out of ' num2str(length(analysed_leaks)) ' ########']);
    hnom = hnom_M(:,index); hleak = hleak_M(:,index);
    qleak = qleak_M(:,index);
    dleak = dleak_M(:,index);
    for mode = 1:2 % 1 - nominal | 2 - leak
        if mode == 1
            hN = hnom(sensors);
        else
            hN = hleak(sensors);
        end

        if strcmp(networkID,'LTOWN')
            hN(sensors==reservoirsID(1)) = 76; % to avoid heads higher than the "reservoirs"
            hN(sensors==reservoirsID(2)) = 76;
        end
        
        % GSI %
        tic
        [x,optinfo_GSI] = GSI(hN,Ld,S,I,reservoirsID,sensors,tau);
        elapsed_time_GSI(i) = toc;

        if mode == 1
            disp('### Nominal scenario started ###');
            x_nom = x;
            x_nom_M(:,index) = x_nom;
            error_nom_GSI_M(:,index) = hnom-x_nom;
            rmse_nom_GSI_M(index) = rmse(hnom,x_nom);
            disp('### Nominal scenario ended ###');
        else
            disp('### Leak scenario started ###');
            x_leak = x;
            x_leak_GSI_M(:,index) = x_leak;
            error_leak_GSI_M(:,index) = hleak-x_leak;
            rmse_leak_GSI_M(index) = rmse(hleak,x_leak);
            [qGSI,~,~] = compute_qandB(x_leak_GSI_M(:,index),E,T);
            rmse_q_GSI(index) = rmse(qGSI,abs(qleak));
            disp('## GSI done ##');
            % AW-GSI %
            [G_aw,B_aw,I2] = compute_GandB(x_nom,d,A,reservoirsID,E,pipes_list);
            [WLap2,Omega,Phi] = compute_LapfromGandB(x_nom,G_aw,B_aw,2,WA,E);
            Phi1 = diag(1./diag(Phi));
            tic
            res_AW = GSI_res(hN-x_nom(sensors),WLap2,S,-I2,reservoirsID,sensors,tau);
            elapsed_time_AWGSI(i) = toc;
            x_leak_AW = res_AW + x_nom;
            x_leak_AWGSI_M(:,index) = x_leak_AW;
            error_leak_AWGSI_M(:,index) = hleak-x_leak_AW;
            rmse_leak_AWGSI_M(index) = rmse(hleak,x_leak_AW);
            [qAWGSI,~,~] = compute_qandB(x_leak_AWGSI_M(:,index),E,T);
            rmse_q_AWGSI(index) = rmse(qAWGSI,abs(qleak));
            disp('## AW-GSI done ##');

            % UKF-AW-GSI %
            ymeas = [hleak(sensors);dleak(idk)];
            R = diag([0.0001*ones(length(sensors),1);0.0001*ones(length(idk),1)]);
            Q = 1*eye(N);
            tic
            [x_leak_UKFAWGSI] = HeadUKF(x_leak_AW,ymeas,R,Q,idk,T,invT,Ss,E,Omega,Phi1);
            elapsed_time_UKFAWGSI(i) = toc;
            x_leak_UKFAWGSI_M(:,index) = x_leak_UKFAWGSI;
            error_leak_UKFAWGSI_M(:,index) = hleak-x_leak_UKFAWGSI;
            rmse_leak_UKFAWGSI_M(index) = rmse(hleak,x_leak_UKFAWGSI);
            [qUKFAWGSI,~,~] = compute_qandB(x_leak_UKFAWGSI_M(:,index),E,T);
            rmse_q_UKFAWGSI(index) = rmse(qUKFAWGSI,abs(qleak));
            disp('## UKF-AW-GSI done ##');

            % D-UKF-AW-GSI %
            known_flows = find_reservoirpipes(E,reservoirsID);
            ymeas_h = [hleak(sensors);dleak(idk)]; %this is completed inside the function
            ymeas_q = [abs(qleak(known_flows))];    %this is completed inside the function
            Rh = diag([0.0001*ones(length(sensors),1);0.0001*ones(length(idk),1); ...
                1000*ones(length(qleak),1)]);
            Rq = diag([1e-6*ones(length(known_flows),1);1e-5*ones(length(qleak),1)]);
            Qh = 1*eye(length(hleak));
            Qq = 1e-5*eye(length(qleak));
            tic
            [x_leak_DUKFAWGSI,q_leak_DUKFAWGSI] = HeadFlowUKF(x_leak_AW,ymeas_h,ymeas_q,known_flows,Rh,Qh,Rq,Qq,idk,T,invT,Ss,E,Omega,Phi1,1);
            elapsed_time_DUKFAWGSI(i) = toc;
            x_leak_DUKFAWGSI_M(:,index) = x_leak_DUKFAWGSI;
            error_leak_DUKFAWGSI_M(:,index) = hleak-x_leak_DUKFAWGSI;
            rmse_leak_DUKFAWGSI_M(index) = rmse(hleak,x_leak_DUKFAWGSI);
            rmse_q_DUKFAWGSI(index) = rmse(abs(qleak),q_leak_DUKFAWGSI);
            disp('## D-UKF-AW-GSI done ##');
            disp('### Leak scenario ended ###');
        end
    end
end

%% Plotting %%

%save Modena_estimation_results.mat rmse_leak_GSI_M rmse_leak_AWGSI_M rmse_leak_UKFAWGSI_M rmse_leak_DUKFAWGSI_M rmse_q_GSI rmse_q_AWGSI rmse_q_UKFAWGSI rmse_q_DUKFAWGSI

% Head %

if enable_plots(cnt_plots)==1
    figure;b = bar([rmse_leak_GSI_M' rmse_leak_AWGSI_M' rmse_leak_UKFAWGSI_M' rmse_leak_DUKFAWGSI_M']);
    b(1).FaceColor = 'b'; % Color for the first group
    b(2).FaceColor = 'r'; % Color for the second group
    b(3).FaceColor = 'k'; % Color for the third group
    b(4).FaceColor = 'g'; % Color for the fourth group
    xlabel('Leak ID','interpreter','latex','fontsize',14);
    ylabel('$RMSE(h,h_{est})$','interpreter','latex','fontsize',14);
    hold on
    hline(mean(rmse_leak_GSI_M),'b--')
    hline(mean(rmse_leak_AWGSI_M),'r--');
    hline(mean(rmse_leak_UKFAWGSI_M),'k--');
    hline(mean(rmse_leak_DUKFAWGSI_M),'g--');legend('GSI','AW-GSI','UKF-AW-GSI','D-UKF-AW-GSI','interpreter','latex','fontsize',10,'NumColumns',2)
end
cnt_plots = cnt_plots + 1;

% Flow %

if enable_plots(cnt_plots)==1
    figure;b = bar([rmse_q_GSI' rmse_q_AWGSI' rmse_q_UKFAWGSI' rmse_q_DUKFAWGSI']);
    b(1).FaceColor = 'b'; % Color for the first group
    b(2).FaceColor = 'r'; % Color for the second group
    b(3).FaceColor = 'k'; % Color for the third group
    b(4).FaceColor = 'g'; % Color for the fourth group
    xlabel('Leak ID','interpreter','latex','fontsize',14);
    ylabel('$RMSE(q,q_{est})$','interpreter','latex','fontsize',14);
    hold on
    hline(mean(rmse_q_GSI),'b--')
    hline(mean(rmse_q_AWGSI),'r--');
    hline(mean(rmse_q_UKFAWGSI),'k--');
    hline(mean(rmse_q_DUKFAWGSI),'g--');legend('GSI','AW-GSI','UKF-AW-GSI','D-UKF-AW-GSI','interpreter','latex','fontsize',10,'NumColumns',2)
end
cnt_plots = cnt_plots + 1;