%% ########### Evaluation of estimation results for UKF-(AW)GSI and D-UKF-AW-GSI at L-TOWN network ########## %%

clear; clc; rng('default')
addpath('functions');

%% Configuration %%

networkID = 'LTOWN';

% Load datasets %

hyd_data_file = '2018_SCADA_Measurements_zoneA_withdemands.mat';
struc_data_file = {'coordenades_BattleDIM1_1.mat'};

enable_plots = [1 1 1]; cnt_plots = 1;
extract_nominal_data = 1;

load(fullfile('data',networkID,'hydraulic',hyd_data_file));
for i=1:length(struc_data_file)
    load(fullfile('data',networkID,'structure',struc_data_file{i}));
end

% Structural data %

N = length(A);

reservoirsID = [111 300];
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

zone_to_study = zoneA;

% Extract pipes of Zone %

E_comp = d.NodesConnectingLinksIndex;

pipes_in_zone = []; Ec = []; k=1;
for pipe=1:length(E_comp)
    if ~isempty(intersect(zone_to_study,E_comp(pipe,1))) && ~isempty(intersect(zone_to_study,E_comp(pipe,2)))
        pipes_in_zone = [pipes_in_zone;pipe];
        E(k,:) = [find(zone_to_study==E_comp(pipe,1)) find(zone_to_study==E_comp(pipe,2))];
        k=k+1;
    end
end

nEdges = length(E);

Gorig = G; Dorig = D; PipeDistanceorig = PipeDistance;

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

sensors_comp = [sensors;reservoirsID_A]; 
sensors_in_zone = intersect(zone_to_study,sensors_comp);
sensors_in_zone_nR = intersect(zone_to_study,sensors);
sensors_nR = zeros(length(sensors_in_zone_nR),1);
sensors_a = zeros(length(sensors_in_zone),1);
for s=1:length(sensors_in_zone)
    sensors_a(s) = find(zone_to_study==sensors_in_zone(s));
end
for s=1:length(sensors_in_zone_nR)
    sensors_nR(s) = find(zone_to_study==sensors_in_zone_nR(s));
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

pressureRaw_orig = pressureRaw;
if extract_nominal_data == 1
    pressureRaw = deletePumpEffect(pressureRaw,flowsRaw);
end


% (optional) Plot network %

% if enable_plots(cnt_plots)==1
%     figure(cnt_plots);
%     fig = plot(G,'XData',node_coordenades(:,1),'YData',node_coordenades(:,2),'Marker','o','MarkerSize',4);
%     xlabel('Latitude','interpreter','latex','fontsize',14);
%     ylabel('Longitude','interpreter','latex','fontsize',14);
%     % title(sprintf('%s - Network graph',network_name));
%     grid
% 
%     highlight(fig,reservoirsID,'MarkerSize',7,'NodeColor','g','Marker','s');
%     highlight(fig,sensors,'MarkerSize',7,'NodeColor','m');%,'Marker','s');
% end
% cnt_plots = cnt_plots + 1;

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
lengths(end-3:end) = 1; roughness(end-3:end) = 140;
pipes_list = pipes_in_zone;

t = 10.674*lengths./(diameters.^4.871.*roughness.^1.852);
T = diag(t(pipes_list)); invT = inv(T);

%% Localization %%

leak_index_v = {461,232,427,628,538,866,654,810,183,158,369};

% Add the 76 m of head to the new reservoirs %
del=[1 2 3 9]; % sensor IDs for the nodes not in Area A
alldel = setdiff(1:size(pressureRaw,2),del);
pressureRaw = pressureRaw(:,alldel)';
hN_total = pressureRaw+d.NodeElevations(sensors_in_zone_nR)';
hN_total = [hN_total(1:2,:);76*ones(1,size(hN_total,2));
    hN_total(3:9,:);76*ones(1,size(hN_total,2));
    hN_total(10:end,:)];

del=[66 230];
alld = 1:659; alldel = setdiff(alld,del);
demandsRaw = demandsRaw'/(1000*3600);

str_nom = 'zoneA_rn_76.mat';
vsensors_file = 'virtual_sensors_100(657nodes).mat';
load(fullfile('data',networkID,'structure',vsensors_file));
idk = vsensors;

st = 6; % sampling step

for leak_index = 1:length(leak_index_v)
    leaky_link_ID = leak_index_v{leak_index};
    switch leaky_link_ID
        case 461
            initial_datetime_v = datetime(2018,3,27,20,35,0);
        case 232
            initial_datetime_v = datetime(2018,2,3,16,5,0);
        case 427
            initial_datetime_v = datetime(2018,6,23,12,0,0);
        case 628
            initial_datetime_v = datetime(2018,5,16,8,0,0);
        case 538
            initial_datetime_v = datetime(2018,5,29,22,0,0);
        case 866
            initial_datetime_v = datetime(2018,6,2,12,0,0);
        case 654
            initial_datetime_v = datetime(2018,9,16,21,5,0);
        case 810
            initial_datetime_v = datetime(2018,11,8,21,25,0);
        case 183
            initial_datetime_v = datetime(2018,8,12,18,30,0);
        case 158
            initial_datetime_v = datetime(2018,10,6,2,40,0);
        case 369
            initial_datetime_v = datetime(2018,11,7,10,25,0);
        otherwise
            error('This leak is not configured.')
    end

for ind = 1:length(leaky_link_ID)
    leaky_link = d.NodesConnectingLinksIndex(leaky_link_ID(ind),:);
    leak(ind,:) = [find(zone_to_study==leaky_link(1)) 
                   find(zone_to_study==leaky_link(2))];
end

% Studied time instants %

initial_datetime = initial_datetime_v;
fin_datetime   =   initial_datetime + days(1);
dtini = abs((datenum(timeRaw)-datenum(initial_datetime))); [~,iniT] = min(dtini);
dtfin = abs((datenum(timeRaw)-datenum(fin_datetime)));     [~,endT] = min(dtfin);

iniT_prevWeek = 1; % index of the time instants considered for the nominal data
endT_prevWeek = iniT_prevWeek+2015;

% Load nominal data / Generate nominal data if necessary %

if extract_nominal_data == 0
    [x_nom_comp,optInfo_nom] = generate_nominal_data(hN_total(:,iniT_prevWeek:endT_prevWeek),N,Ld,S,I,reservoirsID,sensors);
    save(fullfile('data',networkID,'hydraulic',['nominal_heads_iniT_prevWeek=' num2str(iniT_prevWeek)...
          '_endT_prevWeek=' num2str(endT_prevWeek) '_' str_nom]),'x_nom_comp','optInfo_nom',"iniT_prevWeek","endT_prevWeek");
    
else
    load(fullfile('data',networkID,'hydraulic',['nominal_heads_iniT_prevWeek=' num2str(iniT_prevWeek)...
          '_endT_prevWeek=' num2str(endT_prevWeek) '_' str_nom]));
end

wday_ini = weekday(initial_datetime);
if wday_ini == 1
    wday_ini = 8;
end
check_shift = circshift(1:2016,2016-(wday_ini-2)*288);

% End the configuration of the nominal data
x_nom_comp = x_nom_comp(:,check_shift);
% demand_nom = demandsRaw(iniT_prevWeek:endT_prevWeek,alldel)';

% Localization loop %
for mode=1:2 % 1 - AW-GSI | 2 - UKF-AW-GSI
Dsum=0; Dsum2=0; Dsum3 = 0;
% textprogressbar('GSI progress: ');
for t = iniT:st:endT
    
    fprintf('\n######### Time instant: %d out of %d #########\n',t-iniT+1,endT-iniT+1);

    % Extract nominal heads time index that corresponds to the same
    % day/time of the leaky data
    
    if mod(t-(iniT-1),2016)==0 % 2016 => amount of #5min in a week
        index_nominal = 2016;
    else
        index_nominal = mod(t-(iniT-1),2016);
    end

    x_nom = x_nom_comp(:,index_nominal);
    x_nom_time(:,t) = x_nom;

    if mode == 1
        % GSI running

        [x_leak,optInfo_leak(t)] = GSI(hN_total(:,t),Ld,S,I,reservoirsID,sensors,1e3);
        x_leak_time(:,t) = x_leak;

        [G_aw,B_aw,I2] = compute_GandB(x_nom,d,A,reservoirsID,E,pipes_in_zone);
        [WLap2,Omega,Phi] = compute_LapfromGandB(x_nom,G_aw,B_aw,2,WA,E);
        Phi1 = diag(1./diag(Phi));
        Phi1_orig = Phi1;
        Omega_orig = Omega;

        [fnomAW,optinfo_fnomAW] = GSI(x_nom(sensors),WLap2,S,-I2,...
            reservoirsID,...
            sensors,tau);

        res_AW = GSI_res(hN_total(:,t)-x_nom(sensors),WLap2,S,-I2,...
            reservoirsID,...
            sensors,tau);
        x_leak = res_AW + x_nom;
        x_leak_time(:,t) = x_leak;

    else

        [x_leak,optInfo_leak(t)] = GSI(hN_total(:,t),Ld,S,I,reservoirsID,sensors,1e3);
        x_leak_time(:,t) = x_leak;

        [G_aw,B_aw,I2] = compute_GandB(x_nom,d,A,reservoirsID,E,pipes_in_zone);
        [WLap2,Omega,Phi] = compute_LapfromGandB(x_nom,G_aw,B_aw,2,WA,E);
        Phi1 = diag(1./diag(Phi));
        Phi1_orig = Phi1;
        Omega_orig = Omega;

        [fnomAW,optinfo_fnomAW] = GSI(x_nom(sensors),WLap2,S,-I2,...
            reservoirsID,...
            sensors,tau);

        res_AW = GSI_res(hN_total(:,t)-x_nom(sensors),WLap2,S,-I2,...
            reservoirsID,...
            sensors,tau);
        x_leak = res_AW + x_nom;
        x_leak_time(:,t) = x_leak;

        invT = inv(T);

        for mode_2 = 1:2

            if mode_2 == 1
                % Your initial state guess at time k, utilizing measurements up to time k-1: xhat[k|k-1]
                initialStateGuess = [x_nom]; 
                demand = demandsRaw(alldel,index_nominal);
                ymeas = [x_nom(sensors_a);demand(idk)];
            else
                % Your initial state guess at time k, utilizing measurements up to time k-1: xhat[k|k-1]
                initialStateGuess = [x_leak];
                demand = demandsRaw(alldel,t);
                ymeas = [hN_total(:,t);demand(idk)];
            end

            R = diag([0.0001*ones(length(sensors),1); ...
                0.0001*ones(length(idk),1)]);

            Q = 1*eye(length(initialStateGuess));
            gamma = length(idk)/length(demand); 

            x = HeadUKF(initialStateGuess,ymeas,R,Q,idk,T,invT,Ss,E,Omega,Phi1,75);

            if mode_2 == 1
                x_nom = x;
            else
                x_leak = x;
            end
        end
    end

    % LCSM running

    festC=fit(x_nom,x_leak,'poly1');
    coeff = coeffvalues(festC);
    q = [coeff(1);-1;coeff(2)]./sqrt(coeff(1)^2 + (-1)^2);
    Pa = [x_nom x_leak ones(length(x_nom),1)];
    
    D = Pa*q; 
    D=D/max(D); 
    Dsum = Dsum + D;

end

%% Results analysis %%

Da = Dsum;
mu_da = mean(Da); ro_da = std(Da);

candidate_nodes = find(Da>mu_da+ro_da);
prob_cands = Da(Da>mu_da+ro_da);
prob_n = prob_cands - min(prob_cands);
prob_cands = prob_n./max(prob_n);

max_cand = cargmax(prob_cands,5);
best_cand = cargmax(prob_cands);

figure;
fig = plot(G,'XData',node_coordenades(:,1),'YData',node_coordenades(:,2),'Marker','o','MarkerSize',4);
xlabel('X Coordinate [m]','interpreter','latex','fontsize',14);
ylabel('Y Coordinate [m]','interpreter','latex','fontsize',14);

highlight(fig,candidate_nodes,'NodeColor','c','MarkerSize',5);
highlight(fig,sensors,'MarkerSize',5,'NodeColor','r');%,'Marker','s');
highlight(fig,reservoirsID,'MarkerSize',7,'NodeColor','g','Marker','s');
hold on;

highlight(fig,sensors(1:end-2),'MarkerSize',5,'NodeColor','r');

hold on;
ax(1) = plot(NaN,NaN,'gs','MarkerFaceColor','g'); label{1} = sprintf('Reservoir');
ax(2) = plot(NaN,NaN,'ro','MarkerFaceColor','r','LineWidth',1); label{2} = sprintf('Sensors');
hold off;
legend(ax,label,'Location','best','FontSize',14);

leak_name=[];
for ind=1:length(leaky_link_ID)
    if ind > 1
        leak_name = [leak_name '_'];
    end
    larray = num2str(leak_index_v{leak_index});
    leak_name = [leak_name larray(ind,:)];
end

G1 = subgraph(G,candidate_nodes);

figure;%(cnt_plots);
fig = plot(G1,'XData',node_coordenades(candidate_nodes,1),'YData',node_coordenades(candidate_nodes,2),'Marker','o','MarkerSize',4);
xlabel('Latitude','interpreter','latex','fontsize',14);
ylabel('Longitude','interpreter','latex','fontsize',14);
fig.NodeLabel={};

colour_down = [0, 0, 0];
colour_up = [0, 1, 1];
vin = prob_cands;
malla = min(vin):0.001:max(vin);
len = length(malla);
% Determine the cutoff for the top 5%
cutoff_index = round(len * 0.3);

% Create color gradient
colors_malla = [linspace(colour_down(1), colour_up(1), len)', ...
    linspace(colour_down(2), colour_up(2), len)', ...
    linspace(colour_down(3), colour_up(3), len)'];

% Assign black color to the lower 95%
colors_malla(1:cutoff_index, :) = repmat(colour_down, cutoff_index, 1);

% Interpolate the top 5% to go from black to yellow
top_colors = [linspace(0, 1, len - cutoff_index)', ...
    linspace(0, 1, len - cutoff_index)', ...
    linspace(0, 0, len - cutoff_index)'];

colors_malla(cutoff_index+1:end, :) = top_colors;
for i=1:length(prob_cands)
    dif = prob_cands(i) - malla;
    ind = cargmin(abs(dif));
    highlight(fig,i,'NodeColor',colors_malla(ind,:),'MarkerSize',4);
end

colormap(colors_malla);
colorbar;
leak_name=[];
for ind=1:length(leaky_link_ID)
    if ind > 1
        leak_name = [leak_name '_'];
    end
    larray = num2str(leak_index_v{leak_index});
    leak_name = [leak_name larray(ind,:)];
end

cnt_plots = cnt_plots + 1;

%% Compute metrics
% for i=1:length(struc_data_file)
%     G = load(fullfile('data',networkID,'structure',struc_data_file{i})).G;
% end
th = 0.7;
vin=Da;
if min(vin) < 0
    vin_m = vin - min(vin);
else
    vin_m = vin + min(vin);
end
vin_n = vin_m/max(vin_m);

zero_nodes = [reservoirsID';66;200];
if mode > 3 % mode > 3 (node-level metrics) | mode < 3 (area-level metrics)
    for j=1:length(E)
        if ~isempty(find(candidate_nodes==E(j,1))) && ~isempty(find(candidate_nodes==E(j,2)))
            vin_pipe(j) = mean(vin_n(E(j,:)));
        else
            vin_pipe(j) = 0;
        end
    end
    candidates = find(vin_pipe>th);
    if isempty(candidates)
        disp('d')
        candidates = find(vin_pipe>0.45);%mean(vin_pipe)+2*std(vin_pipe));
    end
else
    for j=1:length(E)
        if isempty(find(zero_nodes==E(j,1))) && isempty(find(zero_nodes==E(j,2)))
            vin_pipe(j) = mean(vin_n(E(j,:)));
        else
            vin_pipe(j) = 0;
        end
    end
    candidates = cargmax(vin_pipe);
end

node_candidates = E(candidates,:);
if numel(candidates)>=2
    vin_exp = vin_pipe(candidates)-min(vin_pipe(candidates));
    vin_exp = vin_exp/max(vin_exp);
    vin_norm = vin_exp/sum(vin_exp);
else
    vin_norm = 1;
end

leak_in_candidates = 0;
% zone_to_study = zoneA;
clearvars meters_to_leak pipes_to_leak 

Gzone = subgraph(Gorig,zone_to_study); Gcands = subgraph(Gzone,candidate_nodes);
Dzone=Dorig(zone_to_study,zone_to_study); PDzone=PipeDistanceorig(zone_to_study,zone_to_study);

% Average candidate-to-leak distance (meters)
general_dist = Gzone.distances;

for j=1:length(candidates)
    meters_to_leak(j) = sum(sum(general_dist(leak,node_candidates(j,:))))/4;
    pipes_to_leak(j) = sum(sum(Dzone(leak,node_candidates(j,:))))/4;
    if isequal(leak,node_candidates(j,:))
        meters_to_leak(j) = 0;
        pipes_to_leak(j) = 0;
        leak_in_candidates = 1;
    end
end

leak_number(leak_index,mode) = leaky_link_ID;
score_meters(leak_index,mode) = meters_to_leak*vin_norm';
score_pipes(leak_index,mode) = pipes_to_leak*vin_norm';
search_area_incands(leak_index,mode) = (length(candidates)/Gcands.numedges)*100;
search_area_inzone(leak_index,mode) = (length(candidates)/Gzone.numedges)*100;
success(leak_index,mode) = leak_in_candidates;
end
%%
stats = [leak_number(leak_index,:); success(leak_index,:); score_meters(leak_index,:); score_pipes(leak_index,:); search_area_incands(leak_index,:); search_area_inzone(leak_index,:)]
end