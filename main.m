%% PMBM and PMB implementations with arbitrary clutter

close all
clear
clc
dbstop if error

rng (2022)

plot_enable = true;

%% Parameters of multi-object dynamic model

%scenario: targets are born from a Gaussian density with covariance covering the entire field-of-view and move randomly.

%survival probability
Ps = 0.99;

%Target kinematic state [x-position,y-position,x-velocity,y-velocity]
%Target extent state [orientation,semi-axis length 1,semi-axis length 2]

%Parameters of a nearly constant velocity motion model
%kinematic state dimension
dxr = 4;
%time interval
Ts = 1;
%transition matrix for kinematic state
Ar = [1 0 Ts 0;
    0 1 0 Ts;
    0 0 1 0;
    0 0 0 1];
%process noise
q = 0.01;
%process noise covariance matrix for kinematic state
Cwr = q*[Ts^3/3 0      Ts^2/2 0;
    0      Ts^3/3 0      Ts^2/2;
    Ts^2/2 0      Ts     0;
    0      Ts^2/2 0      Ts];
%measurement rate parameter used for prediction of gamma distribution
eta = 1.05;
%forgetting factor used for prediction of inverse-Wishart distribution
tau = 20;

%struct representation
motionmodel.Ps = Ps;
motionmodel.Ts = Ts;
motionmodel.dxr = dxr;
motionmodel.Ar = Ar;
motionmodel.Cwr = Cwr;
motionmodel.eta = eta;
motionmodel.tau = tau;

%Poisson birth model with Gaussian intensity
%Poisson intensity
birthmodel = struct('w',log(0.1),'xr',[],'Cr',diag([50 50 1 1])^2,'V',[],'v',100,'alpha',500,'beta',100);
%specify kinematic state (x-position,y-position,x-velocity,y_velocity)
birthmodel.xr = [150 150 0 0]';
%specify extent state (orientation,two axis lengths)
birthmodel.V = (100+6)*diag([4 4]);
%Poisson birth rate
b_lambda_0 = 5;
b_lambda = exp(birthmodel.w);

%% Parameters of multi-object measurement model

%detection probability
Pd = 0.9;
%Poisson clutter rate
c_lambda = 20;

%surveillance area
range_x = [0 300];
range_y = [0 300];

%PPP clutter with uniform density
c_pdf = 1/(range_x(2)-range_x(1))/(range_y(2)-range_y(1));
%Poisson clutter intensity
c_intensity = c_lambda*c_pdf;

%stationary clutter source
%number of sources
n_s_clutter = 4;
stationary_clutter = repmat(struct('p_d_clutter',0.98,'lambda_clutter',10,'x_clutter',[],'X_clutter',diag([2 2])),n_s_clutter,1);
stationary_clutter(1).x_clutter = [100;100];
stationary_clutter(2).x_clutter = [100;200];
stationary_clutter(3).x_clutter = [200;100];
stationary_clutter(4).x_clutter = [200;200];

%Parameters of a linear Gaussian measurement model
%measurement dimension
dz = 2;
%observation matrix
H = [1 0 0 0;0 1 0 0];
%covariance of the multiplicative noise
Ch = diag([1/4 1/4]);
%covariance of the measurement noise
Cv = diag([1/4 1/4]);

%struct representation
measmodel.dz = dz;
measmodel.H = H;
measmodel.Ch = Ch;
measmodel.Cv= Cv;
measmodel.Pd = Pd;
measmodel.c_intensity = c_intensity;

%% Generating ground truth

%total time steps
T = 81;

%create memory to store ground truth
%birth time, death time, target state
gt = struct('x_bt',[],'x_dt',[],'x_lambda',[],'xr',[],'X',[]);

%number of targets
nt = 0;
for i = 1:T-1
    %sample the number of newborn targets
    nb = poissrnd(b_lambda);
    %for the first time step, make sure that at least one target is born
    if i == 1
        nb = poissrnd(b_lambda_0);
        if nb == 0
            nb = 1;
        end
    end
    for j = 1:nb
        %number of targets increases by 1
        nt = nt + 1;
        %sample a Gaussian component in the Poisson birth intensity
        b_idx = find(rand < cumsum([0 exp([birthmodel.w])]/b_lambda),1) - 1;
        %sample an initial target state
        gt(nt).x_bt = i;
        gt(nt).x_dt = i;
        %assume fixed Poisson rate
        gt(nt).x_lambda = gamrnd(birthmodel(b_idx).alpha,1/birthmodel(b_idx).beta);
        gt(nt).xr = mvnrnd(birthmodel(b_idx).xr',birthmodel(b_idx).Cr)';
        gt(nt).X = iwishrnd(birthmodel(b_idx).V,birthmodel(b_idx).v-3);
    end
    %generate target trajectory for all the newborn targets
    for j = 1:nt
        %termintate the trajectory if the target dies or moves out of the surveillance area
        %also assumes that no target dies if there is only one
        if (rand < Ps || nt==1) && gt(j).xr(1,end) >= range_x(1) && gt(j).xr(1,end) <= range_x(2) && gt(j).xr(2,end) >= range_y(1) && gt(j).xr(2,end) <= range_y(2) && gt(j).x_dt == i
            gt(j).x_dt = i+1;
            gt(j).x_lambda = [gt(j).x_lambda gt(j).x_lambda(end)];
            %add motion noise when generating trajectory
            gt(j).xr = [gt(j).xr mvnrnd((Ar*gt(j).xr(:,end))',Cwr)'];
            gt(j).X = cat(3,gt(j).X,gt(j).X(:,:,end));
        end
    end
end

%cardinality of multi-target states
card = zeros(T,1);
for i = 1:T
    for j = 1:nt
        if gt(j).x_bt <= i && gt(j).x_dt >= i
            card(i) = card(i) + 1;
        end
    end
end

%% Plot ground truth

if plot_enable
    figure
    grid on
    hold on
    for i = 1:nt
        %plot trajectory
        plot(gt(i).xr(1,1),gt(i).xr(2,1),'b','Marker','o','MarkerFaceColor','b');
        text(gt(i).xr(1,1)+5,gt(i).xr(2,1),int2str(gt(i).x_bt))
        text(gt(i).xr(1,end)-10,gt(i).xr(2,end),int2str(gt(i).x_dt),'Color','red')
        gt_ks_plot = plot(gt(i).xr(1,:),gt(i).xr(2,:),'b');
        %plot extent for every 10 time steps
        for j = 1:10:(gt(i).x_dt - gt(i).x_bt + 1)
            gt_es_plot = plot_extent_iw(gt(i).xr(1:2,j),gt(i).X(:,:,j),'-','b',1);
        end
    end
    %plot stationary source
    for i = 1:n_s_clutter
        sc_plot = plot_extent_iw(stationary_clutter(i).x_clutter,stationary_clutter(i).X_clutter,'-','k',2);
    end
    xlim(range_x)
    ylim(range_y)
    xlabel('x position (m)')
    ylabel('y position (m)')
%     legend([gt_ks_plot, gt_es_plot, sc_plot], {'Target trajectory','Target extent','Stationary source'},'location','best');
end

%% Generate measurements

Z = cell(T,1);
%target-generated measurements
for i = 1:nt
    for t = gt(i).x_bt:gt(i).x_dt
        %generate measurements only if the target is detected
        %no misdetection at time 1
        if rand < Pd || t == 1
            nz = poissrnd(gt(i).x_lambda(t-gt(i).x_bt+1));
            Z{t} = [Z{t} mvnrnd(gt(i).xr(1:2,t-gt(i).x_bt+1)',Ch*gt(i).X(:,:,t-gt(i).x_bt+1)+Cv,nz)'];
        end
    end
end

%append Poisson clutter
for i = 1:T
    nz = poissrnd(c_lambda);
    zx = rand(nz,1)*(range_x(2)-range_x(1)) + range_x(1);
    zy = rand(nz,1)*(range_y(2)-range_y(1)) + range_y(1);
    Z{i} = [Z{i} [zx zy]'];
end

%append stationary clutter
for i = 1:T
    for j = 1:n_s_clutter
        if rand < stationary_clutter(j).p_d_clutter
            nz = poissrnd(stationary_clutter(j).lambda_clutter);
            zc = mvnrnd(stationary_clutter(j).x_clutter',stationary_clutter(j).X_clutter,nz)';
            Z{i} = [Z{i} zc];
        end
    end
end

%% Run filter

%paramter setting

%gate size in probability
% paras.gating.Pg = 0.999;
% paras.gating.size = chi2inv(paras.gating.Pg,dz);
paras.gating.size = 20;
%hyperparameters in DBSCAN
paras.dbscan.max_dist = 5;
paras.dbscan.min_dist = 0.2;
paras.dbscan.grid_dist = 0.1;

%pruning threshold for global hypotheses
paras.pruning.w = log(1e-2);
%pruning threshold for ppp intensity
paras.pruning.ppp = log(1e-3);
%merging threshold for ppp intensity
paras.merging.ppp = 2;
%cap of global hypotheses
paras.cap.w = 20;
%pruning threshold for Bernoulli
paras.pruning.r = 1e-3;
%whether to perform recycling
paras.recycle = false;
if paras.recycle
    paras.pruning.r = 1e-1;
end

%whether to perform MB approximation
paras.mb_approx = true;
%M-best assignments in Murty
paras.M = 20;

%estimator to extract multi-target state
paras.estimator = 1;
%estimator 1: extract state from global hypothesis with the highest weight and Bernoulli components with large enough probability of existence
%estimator 2: MAP cardinality estimator
%threshold to extract state estimate from Bernoulli
paras.estimate.r = 0.4;

%parameters of GOSPA metric
gospa.p = 2;
gospa.c = 20;
gospa.alpha = 2;

%initialisation parameters
%global hypothesis weight in logarithm
mbm.w = 0;
%global hypothesis look-up table
mbm.table = zeros(1,0);
%local hypothesis trees (collections of single target hypotheses)
mbm.track = cell(0,1);
%each Bernoulli is parameterised by 1) existence probability r, 2) mean and covariance of the kinematic state xr, Cr, 3) parameters of the extent state V, v, 4) parameters of gamma distribution alpha, beta.
%PPP for undetected targets, initialised using birth model
ppp = birthmodel;

%memory to store state estimate
est = cell(T,1);
card_est = zeros(T,1);
t_elapsed = zeros(T,1);
%recursive Bayesian estimation
fprintf('Time step: ')
for t = 1:T

    fprintf('%d ',t)
    tic

    %ellipsoidal gating for detected targets
    %use a boolean vector to store the gating result of each local hypothesis
    gating_matrix_d = cellfun(@(x) cell2mat(arrayfun(@(x) ellips_gating(x,Z{t},measmodel,paras.gating.size),x,'uniformoutput',false)),mbm.track,'uniformoutput',false);

    %ellipsoidal gating for undetected targets
    %use a boolean vector to store the gating result of each component
    gating_matrix_u = cell2mat(arrayfun(@(x) ellips_gating(x,Z{t},measmodel,paras.gating.size),ppp,'uniformoutput',false)');

    %ellipsoidal gating for stationary sources
    %use a boolean vector to store the gating result of each stationary source
    gating_matrix_s = cell2mat(arrayfun(@(x) ellips_gating_clutter(Z{t},x,paras.gating.size),stationary_clutter,'uniformoutput',false)');

    %remove unused measurements according to the gating result
    %gating result of all targets
    %number of tracks
    n_track = length(mbm.track);
    if n_track > 0
        gating = logical(sum(gating_matrix_u,2) + sum(gating_matrix_s,2) + sum(cell2mat(gating_matrix_d'),2));
    else
        gating = logical(sum(gating_matrix_u,2) + sum(gating_matrix_s,2));
    end

    %used measurements
    W = Z{t}(:,gating);
    %number of measurements after gating
    nm = size(W,2);

    %reconstruct gating matrix
    gating_matrix_d = cellfun(@(x) x(gating,:),gating_matrix_d,'uniformoutput',false);
    gating_matrix_u = gating_matrix_u(gating,:);
    gating_matrix_s = gating_matrix_s(gating,:);
    %use DBSCAN to obtain multiple partitions
    partitions = gen_partitions(W,paras.dbscan,n_track);
    %number of partitions
    np = length(partitions);

    %find all the unique clusters in all measurement partitions
    [clusters,IA,IC] = unique(cell2mat(partitions')','rows');
    clusters = clusters';
    n_clusters = length(IA);
    %reconstruct partitions to let it contain indices of clusters
    nc_p = cellfun(@(x) size(x,2),partitions);
    partitions_indices = cell(np,1);
    idx = 0;
    for i = 1:np
        partitions_indices{i} = IC(idx+1:idx+nc_p(i));
        idx = idx + nc_p(i);
    end

    %create clutter hypotheses for stationary sources
    stationary_hypo = repmat(struct('c',[],'lik',[]),n_s_clutter,1);
    for i = 1:n_s_clutter
        %no clutter hypothesis
        stationary_hypo(i).c = 0;
        stationary_hypo(i).lik = log(1-stationary_clutter(i).p_d_clutter+stationary_clutter(i).p_d_clutter*exp(-stationary_clutter(i).lambda_clutter));
        %create a clutter hypothesis for each clutter
        for c = 1:n_clusters
            %check if the cth clutter is in the gate
            if sum(gating_matrix_s(:,i)-clusters(:,c)<0) == 0
                stationary_hypo(i).c(end+1) = c;
                stationary_hypo(i).lik(end+1) = log(stationary_clutter(i).p_d_clutter)-stationary_clutter(i).lambda_clutter+clutter_lik(W(:,clusters(:,c)),stationary_clutter(i));
            end
        end
    end

    %create single target hypotheses for newly detected targets
    bern_new = repmat(struct('r',0,'xr',zeros(dxr,1),'Cr',ones(dxr,dxr),'V',zeros(2,2),'v',0,'alpha',1,'beta',1),1,n_clusters);
    lik_new = zeros(n_clusters,1);
    for c = 1:n_clusters
        %check if the cth cluster is in the gate of any ppp components
        ppp_idx = sum(gating_matrix_u-clusters(:,c)<0) == 0;
        if any(ppp_idx)
            [lik_new(c),bern_new(c)] = ppp_upd(ppp(ppp_idx),W(:,clusters(:,c)),measmodel);
        else
            %if not, they are all Poisson clutter
            lik_new(c) = sum(clusters(:,c))*log(measmodel.c_intensity);
        end
    end

    %create updated single target hypotheses for detected targets
    %number of single target hypotheses per track
    n_local_hypo = cellfun('length',mbm.track);
    tracks_upd = cell(n_track,1);
    for i = 1:n_track
        tracks_upd{i} = cell(n_local_hypo(i),1);
        for j = 1:n_local_hypo(i)
            %misdetection for the jth single target hypothesis under the ith local hypothesis tree
            [l_missed,bern_missed] = bern_miss(mbm.track{i}(j),measmodel);
            tracks_upd{i}{j}.c = 0;
            tracks_upd{i}{j}.lik = l_missed;
            tracks_upd{i}{j}.bern = bern_missed;
            %measurement update for the jth single target hypothesis under the ith local hypothesis tree
            for c = 1:n_clusters
                %check if the cth cluster is in the gate of the corresponding single target hypothesis
                if sum(gating_matrix_d{i}(:,j)-clusters(:,c)<0) == 0
                    [l_upd,bern_updated] = bern_upd(mbm.track{i}(j),W(:,clusters(:,c)),measmodel);
                    tracks_upd{i}{j}.c(end+1) = c;
                    tracks_upd{i}{j}.lik(end+1) = l_upd;
                    tracks_upd{i}{j}.bern(end+1) = bern_updated;
                end
            end
        end
    end

    tracks_new = cell(n_clusters,1);
    for c = 1:n_clusters
        tracks_new{c}.c = [0 c];
        tracks_new{c}.lik = [0 lik_new(c)];
        tracks_new{c}.bern = struct('r',0,'xr',zeros(dxr,1),'Cr',ones(dxr,dxr),'V',zeros(2,2),'v',0,'alpha',1,'beta',1);
        tracks_new{c}.bern = [tracks_new{c}.bern bern_new(c)];
    end

    %reset m to the number of new tracks
    m = length(tracks_new);

    %update local hypothesis trees
    mbm_upd.track = cell(n_track+m,1);
    for i = 1:n_track
        idx = 0;
        for j = 1:length(tracks_upd{i})
            mbm_upd.track{i} = [mbm_upd.track{i} tracks_upd{i}{j}.bern];
            %use an extra variable to record the index of each new single target hypothesis in local hypothesis tree i
            n_ij = length(tracks_upd{i}{j}.c);
            tracks_upd{i}{j}.idx = (1:n_ij) + idx;
            idx = idx + n_ij;
        end
    end
    for i = 1:m
        mbm_upd.track{n_track+i,1} = tracks_new{i}.bern;
    end

    %data association for each global association hypothesis
    mbm_upd.w = [];
    mbm_upd.table = zeros(0,n_track+m);
    A = length(mbm.w);
    for a = 1:A
        a_indices = mbm.table(a,:);
        costs_temp = 0;
        for j = 1:n_track
            if a_indices(j) > 0
                single_hypo = tracks_upd{j}{a_indices(j)};
                costs_temp = costs_temp + single_hypo.lik(1);
            end
        end
        %if there is no measurement partition, all the targets are misdetected
        if n_clusters==0
            table_upd = zeros(1,n_track+m);
            w_upd = 0;
            for j = 1:n_track
                if a_indices(j) > 0
                    table_upd(1,j) = tracks_upd{j}{a_indices(j)}.idx(1);
                    w_upd = w_upd + tracks_upd{j}{a_indices(j)}.lik(1);
                end
            end
            for j = n_track+1:n_track+m
                table_upd(1,j) = 1;
            end
            mbm_upd.w = [mbm_upd.w;w_upd + mbm.w(a)];
            mbm_upd.table = [mbm_upd.table;table_upd];
        end
        %go through each measurement partition
        for p = 1:np
            %find all the clusters under this partition
            p_idx = partitions_indices{p};
            %number of clusters under this partition
            p_n = length(p_idx);
            %construct cost matrix
            C = inf(p_n,n_track+m);
            for j = 1:n_track
                if a_indices(j) > 0
                    single_hypo = tracks_upd{j}{a_indices(j)};
                    %set cost for detected targets
                    [LIA,LOCB] = ismember(single_hypo.c,p_idx);
                    C(LOCB(LOCB>0),j) = -single_hypo.lik(LIA)' + single_hypo.lik(1);
                end
            end
            %set cost for undetected targets
            for j = n_track+1:n_track+m
                single_hypo = tracks_new{j-n_track};
                [LIA,LOCB] = ismember(single_hypo.c,p_idx);
                C(LOCB(LOCB>0),j) = -single_hypo.lik(LIA)';
            end

            %construct cost matrix for stationary sources
            costs_temp_ac =  costs_temp;
            C_s = inf(p_n,n_s_clutter);
            for j = 1:n_s_clutter
                [LIA,LOCB] = ismember(stationary_hypo(j).c,p_idx);
                if any(LOCB>0)
                    C_s(LOCB(LOCB>0),j) = -stationary_hypo(j).lik(LIA)' + stationary_hypo(j).lik(1);
                    costs_temp_ac = costs_temp_ac + stationary_hypo(j).lik(1);
                end
            end

            %combine cost matrices
            C = [C C_s];

            %find columns of C that contain finite entries
            idx = find(sum(isfinite(C),1) > 0);
            C = C(:,idx);
            %find M-best assignments using Murty's algorithm
            [assignments,~,costs] = kBest2DAssign(C,ceil(paras.M*exp(mbm.w(a))));
            assignments = assignments';
            costs = costs';
            %number of assignments
            n_a = size(assignments,1);
            %restore track indices
            for i = 1:n_a
                assignments(i,:) = idx(assignments(i,:));
            end

            table_upd = zeros(n_a,n_track+m);
            %update the glocal hypothesis look-up table and weight
            for i = 1:n_a
                %go through each association in a given assignment
                for j = 1:p_n
                    %check if detected targets or undetected targets, otherwise cluster is assigned to stationary sources
                    assoc = assignments(i,j);
                    if assoc <= n_track
                        %find the index of the corresponding cluster
                        table_upd(i,assoc) = tracks_upd{assoc}{a_indices(assoc)}.idx(tracks_upd{assoc}{a_indices(assoc)}.c == p_idx(j));
                    elseif assoc <= n_track+m
                        %find the index of the corresponding cluster
                        table_upd(i,assoc) = find(tracks_new{assoc-n_track}.c == p_idx(j),1);
                    end
                end
                %go through each unassociated track
                unassign = true(n_track+m,1);
                unassign(assignments(i,assignments(i,:)<=n_track+m)) = false;
                unassign = find(unassign);
                for j = 1:length(unassign)
                    if unassign(j) <= n_track
                        temp = a_indices(unassign(j));
                        if temp > 0
                            table_upd(i,unassign(j)) = tracks_upd{unassign(j)}{temp}.idx(1);
                        else
                            table_upd(i,unassign(j)) = 0;
                        end
                    else
                        table_upd(i,unassign(j)) = 1;
                    end
                end
            end

            %when computing the global hypothesis weight, assume that each measurement partition is equal likely, i.e., uniform distributed
            mbm_upd.w = [mbm_upd.w;-costs'+costs_temp_ac+mbm.w(a)];
            mbm_upd.table = [mbm_upd.table;table_upd];
        end
    end

    %normalise global hypothesis weight
    mbm_upd.w = normalizeLogWeights(mbm_upd.w);
    %prune updated global hypotheses with small weights
    [mbm_upd_w,order] = sort(exp(mbm_upd.w),'descend');
    pos = find(cumsum(mbm_upd_w)>=1-exp(paras.pruning.w),1);
    mbm_upd.w = mbm_upd.w(order(1:pos));
    mbm_upd.table = mbm_upd.table(order(1:pos),:);
    mbm_upd.w = normalizeLogWeights(mbm_upd.w);
    if ~paras.mb_approx
        %cap the number of global hypotheses
        if length(mbm_upd.w) > paras.cap.w
            [~,idx] = sort(mbm_upd.w,'descend');
            mbm_upd.w = mbm_upd.w(idx(1:paras.cap.w));
            mbm_upd.w = normalizeLogWeights(mbm_upd.w);
            mbm_upd.table = mbm_upd.table(idx(1:paras.cap.w),:);
        end
    end

    %remove single target hypotheses with small probability of existence
    for i = 1:length(mbm_upd.track)
        idx = find(([mbm_upd.track{i}.r] < exp(paras.pruning.ppp)) | ([mbm_upd.track{i}.alpha]./[mbm_upd.track{i}.beta] < 1) | ([mbm_upd.track{i}.v] < 7));
        for j = 1:length(idx)
            mbm_upd.table(mbm_upd.table(:,i) == idx(j),i) = 0;
        end
        %re-index
        idx_0 = mbm_upd.table(:,i) > 0;
        [idx,~,temp] = unique(mbm_upd.table(idx_0,i));
        mbm_upd.table(idx_0,i) = temp;
        mbm_upd.track{i} = mbm_upd.track{i}(idx);
    end
    %remove empty track
    idx = ~cellfun('isempty',mbm_upd.track);
    mbm_upd.track = mbm_upd.track(idx);
    mbm_upd.table = mbm_upd.table(:,idx);
    if isempty(mbm_upd.table)
        mbm_upd.table = zeros(1,0);
        mbm_upd.track = cell(0,1);
        mbm_upd.w = 0;
    end

    %merge rows of global hypothesis look-up table that are the same
    if length(mbm_upd.w) > 1
        [mbm_upd.table,~,IC] = unique(mbm_upd.table,'rows');
        n_a = size(mbm_upd.table,1);
        temp = zeros(n_a,1);
        for i = 1:n_a
            [~,temp(i)] = normalizeLogWeights(mbm_upd.w(IC==i));
        end
        mbm_upd.w = temp;
    end

    %misdetection update of ppp
    ppp = ppp_miss(ppp,measmodel);

    %number of global hypotheses and tracks
    n_mb = size(mbm_upd.table,1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if paras.mb_approx && n_mb > 1
        %find tracks with only one single target hypothesis being included in any of the global hypotheses
        idx = sum(mbm_upd.table - 1,1) ~= 0;
        mbm_upd.track = [mbm_upd.track(idx);mbm_upd.track(~idx)];
        mbm_upd.table = [mbm_upd.table(:,idx) mbm_upd.table(:,~idx)];

        %track-oriented merging
        %number of tracks with conflicts
        n_track_t = length(idx);
        for i = 1:n_track_t
            %compute marginal probability that the single target hypothesis in track i is included in the global hypothesis
            c = unique(mbm_upd.table(:,i));
            c = c(c>0);
            lc = length(c);
            w_margin = zeros(lc,1);
            for j = 1:lc
                %find all the global hypotheses that contain this single target hypothesis and compute the log sum of the weights
                [~,log_sum_w] = normalizeLogWeights(mbm_upd.w(mbm_upd.table(:,i)==c(j)));
                %also take into account the probability of existence
                w_margin(j) = log_sum_w + log(mbm_upd.track{i}(c(j)).r);
            end
            %perform moment matching
            [w_margin_n,log_sum_w] = normalizeLogWeights(w_margin);
            mb_approx(i).r = exp(log_sum_w);
            [mb_approx(i).xr,mb_approx(i).Cr] = kinematic_merge(mbm_upd.track{i},w_margin_n);
            [mb_approx(i).V,mb_approx(i).v] = extent_merge(mbm_upd.track{i},w_margin_n);
            [mb_approx(i).alpha,mb_approx(i).beta] = gamma_merge(mbm_upd.track{i},w_margin_n);
        end

        %track-oriented merging
        for i = 1:n_track_t
            mbm_upd.track{i} = mb_approx(i);
        end

        mbm_upd.w = 0;
        mbm_upd.table = ones(1,n_track_t);

        %recycle tracks with small probability of existence
        idx_logical = (cellfun(@(x) x.r, mbm_upd.track) >= paras.pruning.r) & (cellfun(@(x) x.alpha, mbm_upd.track)./cellfun(@(x) x.beta, mbm_upd.track) > 1) & (cellfun(@(x) x.v, mbm_upd.track) > 7);
        if paras.recycle
            idx = find(~idx_logical);
            for j = 1:length(idx)
                ppp(end+1,1).w = log(mbm_upd.track{idx(j)}.r);
                ppp(end,1).xr = mbm_upd.track{idx(j)}.xr;
                ppp(end,1).Cr = mbm_upd.track{idx(j)}.Cr;
                ppp(end,1).V = mbm_upd.track{idx(j)}.V;
                ppp(end,1).v = mbm_upd.track{idx(j)}.v;
                ppp(end,1).alpha = mbm_upd.track{idx(j)}.alpha;
                ppp(end,1).beta = mbm_upd.track{idx(j)}.beta;
            end
        end
        %delete these single target hypotheses
        mbm_upd.track = mbm_upd.track(idx_logical);
        mbm_upd.table = mbm_upd.table(idx_logical);
        if isempty(mbm_upd.table)
            mbm_upd.table = zeros(1,0);
            mbm_upd.track = cell(0,1);
            mbm_upd.w = 0;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    else
        %recycle single target hypotheses with small probability of existence
        for i = 1:length(mbm_upd.track)
            idx = find(([mbm_upd.track{i}.r] < paras.pruning.r) | ([mbm_upd.track{i}.alpha]./[mbm_upd.track{i}.beta] < 1) | ([mbm_upd.track{i}.v] < 7));
            %for each single target hypothesis to be recycled, find global
            %hypotheses that include it
            for j = 1:length(idx)
                temp = mbm_upd.table(:,i)==idx(j);
                mbm_upd.table(temp,i) = 0;
                if paras.recycle
                    [~,log_sum_w] = normalizeLogWeights(mbm_upd.w(temp));
                    ppp(end+1,1).w = log_sum_w+log(mbm_upd.track{i}(idx(j)).r);
                    ppp(end,1).xr = mbm_upd.track{i}(idx(j)).xr;
                    ppp(end,1).Cr = mbm_upd.track{i}(idx(j)).Cr;
                    ppp(end,1).V = mbm_upd.track{i}(idx(j)).V;
                    ppp(end,1).v = mbm_upd.track{i}(idx(j)).v;
                    ppp(end,1).alpha = mbm_upd.track{i}(idx(j)).alpha;
                    ppp(end,1).beta = mbm_upd.track{i}(idx(j)).beta;
                end
            end
            %delete these single target hypotheses
            %re-index
            idx_0 = mbm_upd.table(:,i) > 0;
            [idx,~,temp] = unique(mbm_upd.table(idx_0,i));
            mbm_upd.table(idx_0,i) = temp;
            mbm_upd.track{i} = mbm_upd.track{i}(idx);
        end
        %remove empty track
        idx = ~cellfun('isempty',mbm_upd.track);
        mbm_upd.track = mbm_upd.track(idx);
        mbm_upd.table = mbm_upd.table(:,idx);
        if isempty(mbm_upd.table)
            mbm_upd.table = zeros(1,0);
            mbm_upd.track = cell(0,1);
            mbm_upd.w = 0;
        end

        %merge rows of global hypothesis look-up table that are same
        if length(mbm_upd.w) > 1
            [mbm_upd.table,~,IC] = unique(mbm_upd.table,'rows');
            n_a = size(mbm_upd.table,1);
            temp = zeros(n_a,1);
            for i = 1:n_a
                [~,temp(i)] = normalizeLogWeights(mbm_upd.w(IC==i));
            end
            mbm_upd.w = temp;
        end
    end

    %remove ppp components with small weights
    idx = ([ppp.w] > paras.pruning.ppp) & ([ppp.alpha]./[ppp.beta] > 1) & ([ppp.v] > 7);
    ppp = ppp(idx);
    if ~paras.mb_approx
        %merge similar ppp components
        ppp = mixtureReduction(ppp,paras.merging.ppp);
    end

    %multi-target state estimation
    if paras.estimator == 1
        %find the global hypothesis with the highest weight
        [~,I] = max(mbm_upd.w);
        est{t}.xr = zeros(dxr,0);
        est{t}.X = zeros(2,2,0);
        if ~isempty(mbm_upd.table)
            hypo_best = mbm_upd.table(I,:);
            for i = 1:length(hypo_best)
                if mbm_upd.table(I,i) > 0
                    %extract state estimate from Bernoulli component with large
                    %enough probability of existence
                    if mbm_upd.track{i}(mbm_upd.table(I,i)).r > paras.estimate.r
                        card_est(t) = card_est(t) + 1;
                        est{t}.xr = [est{t}.xr mbm_upd.track{i}(mbm_upd.table(I,i)).xr];
                        est{t}.X = cat(3,est{t}.X,mbm_upd.track{i}(mbm_upd.table(I,i)).V/(mbm_upd.track{i}(mbm_upd.table(I,i)).v-6));
                    end
                end
            end
        end
        est{t}.card = size(est{t}.xr,2);
    elseif paras.estimator == 2
        %compute the cardinality distribution
        [n_a,n_track] = size(mbm_upd.table);
        card_dist = zeros(1,n_track+1);
        pcard = zeros(n_a,n_track+1);
        for i = 1:n_a
            r = [];
            for j = 1:n_track
                if mbm_upd.table(i,j) > 0
                    r = [r mbm_upd.track{j}(mbm_upd.table(i,j)).r];
                end
            end
            if ~isempty(r)
                %avoid numerical underflow
                r(r>1-1e-6) = 1-1e-6;
                pcard(i,1:length(r)+1) = prod(1-r)*poly(-r./(1-r));
            end
            card_dist = card_dist + pcard(i,:)*exp(mbm_upd.w(i));
        end
        %obtain the maximum cardinality
        [~,card_max] = max(card_dist);
        %find the global hypothesis with the highest weight and the same
        %MAP cardinality estimate
        [~,a_best] = max(pcard(:,card_max));
        r = zeros(n_track,1);
        xr = zeros(dxr,n_track);
        X = zeros(2,2,n_track);
        for i = 1:n_track
            if mbm_upd.table(a_best,i) > 0
                r(i) = mbm_upd.track{i}(mbm_upd.table(a_best,i)).r;
                xr(:,i) = mbm_upd.track{i}(mbm_upd.table(a_best,i)).xr;
                X(:,:,i) = mbm_upd.track{i}(mbm_upd.table(a_best,i)).V/(mbm_upd.track{i}(mbm_upd.table(a_best,i)).v-6);
            end
        end
        [~,I] = sort(r,'descend');
        card_est(t) = card_max-1;
        est{t}.xr = xr(:,I(1:card_max-1));
        est{t}.X = X(:,:,I(1:card_max-1));
        est{t}.card = sum(card_dist.*(0:length(card_dist)-1));
    end

    %prediction step
    %prediction for detected targets
    mbm.w = mbm_upd.w;
    mbm.table = mbm_upd.table;
    mbm.track = cellfun(@(x) arrayfun(@(x) bern_pred(x,motionmodel),x),mbm_upd.track,'uniformoutput',false);

    %prediction for undetected targets
    ppp = ppp_pred(ppp,motionmodel,birthmodel);

    t_elapsed(t) = toc;
end
fprintf('\n')

card_est = cellfun(@(x) x.card,est);
card_err = card_est - card;
%evaluate multi-target filtering performance using GOSPAs
d_gospa = zeros(T,1);
decomposed_cost = repmat(struct('localisation',[],'missed',[],'false',[]),T,1);
for t = 1:T
    x_mat.x = zeros(2,0);
    x_mat.X = zeros(2,2,0);
    for i = 1:length(gt)
        if gt(i).x_bt <= t && gt(i).x_dt >= t
            x_mat.x = [x_mat.x gt(i).xr(1:2,t-gt(i).x_bt+1)];
            x_mat.X = cat(3,x_mat.X,gt(i).X(:,:,t-gt(i).x_bt+1));
        end
    end
    y_mat.x = est{t}.xr;
    y_mat.X = est{t}.X;
    [d_gospa(t),~,decomposed_cost(t)] = GOSPA_extended(x_mat,y_mat,gospa.p,gospa.c,gospa.alpha);
end

if plot_enable

    figure
    plot(1:T,card,'linewidth',2)
    grid on
    hold on
    plot(1:T,card_est,'linewidth',2)
    xlabel('Time step')
    ylabel('Number of targets')
    legend('True cardinality','Estimated cardinality')

    figure
    subplot(2,2,1)
    plot(1:T,sqrt(d_gospa),'linewidth',2)
    grid on
    xlabel('Time step')
    ylabel('RMS GOSPA error')
    subplot(2,2,2)
    plot(1:T,sqrt([decomposed_cost.localisation]),'linewidth',2)
    grid on
    xlabel('Time step')
    ylabel('RMS State estimation error')
    subplot(2,2,3)
    plot(1:T,sqrt([decomposed_cost.missed]),'linewidth',2)
    grid on
    xlabel('Time step')
    ylabel('RMS Missed target error')
    subplot(2,2,4)
    plot(1:T,sqrt([decomposed_cost.false]),'linewidth',2)
    grid on
    xlabel('Time step')
    ylabel('RMS False target error')

end

fprintf('Mean RMS GOSPA: %.2f\n',sqrt(mean(d_gospa)))
fprintf('Total runtime: %.2f\n',sum(t_elapsed))

