clear

num_models = 6; % Number of models to be tested


%%%%%%%%%%  Prepare simulations


% Create a list of colours for consistency
setup_colours = [0.4660 0.6740 0.1880; 0.3010 0.7450 0.9330; 0.9290 0.6940 0.1250; 0.4940 0.1840 0.5560; 1 0 0; 0.5 0.5 0.5];

% Names of models
models = ["Empirical","ER","SBM","Deg-corr SBM","Configuration","Comm-config","Analytic"];
models2 = ["Empirical","ER","SBM/Homogenised","Deg-corr SBM","Configuration","Comm-config"];

num_runs_per_model = 10;  % 1000
maxT = 1500;              % Total number of timesteps, 1500
num_timepoints = 100;     % Number of reccorded timesteps, 250
times = linspace(1,maxT,num_timepoints);


%%%%%%%%%% Get features of empirical network


create_empirical_network
% This provides
%{
- G_empirical : the empirical network
- A_empirical : the empirical adjacency matrix
- deg_empirical : the empirical degree sequence
- deg_w_empirical : the emprical within-degree sequence
- deg_b_empirical : the emprical between-degree sequence
- C0,C1 : community sizes
- deg_corr_SBM_expected_adjacency : expected adjancecy matrix for
  degree-corrected SBM
%}

[m,~,~] = get_network_mixing_rate(G_empirical,C0,C1);
x_analytic = c1*(c1 + c2.*exp(-m(1).*times./N));


%%%%%%%%%%  Prepare arrays for storing outputs


m_store = zeros(1,num_models);

ave_x_store = zeros(num_timepoints,num_models);
ave_node_state = zeros(num_timepoints,N);

observed_mean_state_C0 = zeros(num_timepoints,num_models);
observed_mean_state_C1 = zeros(num_timepoints,num_models);

expected_state_C0 = zeros(num_timepoints,num_models);
expected_state_C1 = zeros(num_timepoints,num_models);

ave_node_state_C0_store = zeros(num_timepoints,C0,num_models);
ave_node_state_C1_store = zeros(num_timepoints,C1,num_models);

expected_adjacency_matrix = zeros(N,N,num_models+1);
deviation_store = zeros(N,num_models);

dynamical_distance_x = zeros(num_models+1,num_models+1);
network_distance = zeros(num_models,num_models);


%%%%%%%%%% Calculate expected adjacency matrix for each model


% 1. Empirical network
expected_adjacency_matrix(:,:,1) = A_empirical;

% 2. ER graph
% Calculate the MLE for the ER model parameter p
p_ER = full(sum(sum(A_empirical))/(N*(N-1)));
% Set the expected adjacency matrix
expected_adjacency_matrix(:,:,2) = p_ER*(1 - eye(N));

% 3. SBM
% Calculate the MLE for SBM parameters
M0 = 0.5*full(sum(sum(A_empirical(1:C0,1:C0))));
M1 = 0.5*full(sum(sum(A_empirical(1+C0:end,1+C0:end))));
Mb = full(sum(sum(A_empirical(1:C0,1+C0:end))));
p_00 = 2*M0/(C0*(C0-1));
p_11 = 2*M1/(C1*(C1-1));
p_01 = Mb/(C0*C1);
% Set the expected adjacency matrix
expected_adjacency_matrix(:,:,3) = p_01;
expected_adjacency_matrix(1:C0,1:C0,3) = p_00*(1 - eye(C0));
expected_adjacency_matrix(1+C0:N,1+C0:N,3) = p_11*(1 - eye(C1));

% 4. Degree-corrected SBM
expected_adjacency_matrix(:,:,4) = deg_corr_SBM_expected_adjacency;

% 5. Configuration model
M_empirical = numedges(G_empirical);
expected_adjacency_matrix(:,:,5) = deg_empirical*deg_empirical'/(2*M_empirical - 1);

% 6. Community configuration model
expected_adjacency_matrix(:,:,6) = [deg_w_empirical(1:C0)'*deg_w_empirical(1:C0)/(2*M0 - 1), deg_b_empirical(1:C0)'*deg_b_empirical(1+C0:N)/(2*Mb - 1);
                                    deg_b_empirical(C0+1:N)'*deg_b_empirical(1:C0)/(2*Mb - 1), deg_w_empirical(1+C0:N)'*deg_w_empirical(1+C0:N)/(2*M1 - 1)];


%%%%%%%%%% Main loop 


for model_setup=1:num_models

    disp(model_setup)
    
    % Store output from each realisation
    states_store = zeros(num_timepoints,N,num_runs_per_model);
    % Create empty array to store the value of m on each run
    m_run = zeros(1,num_runs_per_model);


%%%%%%%%%% Repeatedly simulate dynamics and store output


    parfor run=1:num_runs_per_model
        
        Gr = graph();
        % Generate the appropriate network
        if model_setup==1
                % 1. Empirical network
                Gr = G_empirical;
        else, if model_setup==2
                % 2. ER
                [Gr,~] = create_ER_network(N,p_ER);
        else, if model_setup==3
                % 3. SBM
                [Gr,~] = create_SBM_network([C0,C1],[p_00,p_11],p_01);
        else, if model_setup==4
                % 4. Degree-corrected SBM
                [Gr,~] = create_degree_corrected_SBM(N,deg_corr_SBM_expected_adjacency);
        else, if model_setup==5
                % 5. Configuration model
                [Gr,~] = create_configuration_model_network(deg_empirical);
        else, if model_setup==6
                % 6. Community-configuration model
                [Gr,~] = create_community_configuration_network(deg_w_empirical,deg_b_empirical,C0);
              end
              end
              end
              end
              end
        end

        % Simulate the dynamics
        [~,states] = simulate_exclusion_process([C0,C1],community,Gr,times);
        % Store the output 
        states_store(:,:,run) = states;
        % Get the value of m
        [m_run(run),~,~] = get_network_mixing_rate(Gr,C0,C1);

    end
    

%%%%%%%%%% Average outputs


    m_store(model_setup) = mean(m_run);
    ave_node_state = mean(states_store,3);

    % Separate ave_node_state into C1 nodes and C2 nodes
    ave_node_state_C0 = ave_node_state(:,1:C0);
    ave_node_state_C1 = ave_node_state(:,1+C0:N);
    
    % Calculate the expected state of each node at each timepoint
    probC1isT1 = x_analytic./c1;
    probC1isT2 = 1 - x_analytic./c1;
    probC2isT1 = c1/c2 - x_analytic./c2;
    probC2isT2 = 1 - c1/c2 + x_analytic./c2;
    
    expected_state_C0(:,model_setup) = 0.*probC1isT1 + 1.*probC1isT2;
    expected_state_C1(:,model_setup) = 0.*probC2isT1 + 1.*probC2isT2;

    % Calculate the deviation from the expected state for each node
    deviation = zeros(N,1);
    for i=1:C0
        deviation(i) = norm(ave_node_state_C0(:,i) - expected_state_C0(:,model_setup)')/maxT;
    end
    for i=1:C1
        deviation(i+C0) = norm(ave_node_state_C1(:,i) - expected_state_C1(:,model_setup)')/maxT;
    end

    % Store deviation values for this model
    deviation_store(:,model_setup) = deviation;

    % Calculate the average state of an average node in each community at
    % each timepoint
    observed_mean_state_C0(:,model_setup) = mean(ave_node_state_C0,2);
    observed_mean_state_C1(:,model_setup) = mean(ave_node_state_C1,2);

    % Calculate the average (observed) value of x
    ave_x = c1.*(1 - observed_mean_state_C0(:,model_setup));
    ave_x_store(:,model_setup) = ave_x;

    ave_node_state_C0_store(:,:,model_setup) = ave_node_state_C0;
    ave_node_state_C1_store(:,:,model_setup) = ave_node_state_C1;

end


%%%%%%%%%%  Prepare figures


figure(1)
clf
tl = tiledlayout(2,num_models,"TileIndexing","columnmajor","TileSpacing","compact");

figure(2)
clf
hold on

figure(3)
clf

figure(4)
clf

figure(5)
clf


%%%%%%%%%% Plot resutls


for model_setup=1:num_models

    figure (1)
    xlabel(tl,'Time','FontSize',12)
    
    % Plot the average node state in C0
    nexttile
    plot1a = plot(times,ave_node_state_C0_store(:,:,model_setup),'Color',[0.00 0.45 0.74],'LineWidth',0.1);
    hold on
    plot1c = plot(times,expected_state_C0(:,model_setup),'k--','LineWidth',1.5);
    plot1b = plot(times,observed_mean_state_C0(:,model_setup),'k','LineWidth',1.5);
    ylim([0 1])
    title(models(model_setup),'FontSize',12)
    if model_setup==1
        ylabel('C0 average state','FontSize',12)
    end

    % Plot the average node state in C1
    nexttile
    plot1d = plot(times,ave_node_state_C1_store(:,:,model_setup),'Color',[0.8500 0.3250 0.098],'LineWidth',0.1);
    hold on
    plot(times,expected_state_C1(:,model_setup),'k--','LineWidth',1.5)
    plot(times,observed_mean_state_C1(:,model_setup),'k','LineWidth',1.5)
    ylim([0 1])
    if model_setup==1
        ylabel('C1 average state','FontSize',12)
    end

    % Plot x
    figure (2)
    if model_setup==1
        plot(times,x_analytic,'k--','LineWidth',1.5,'DisplayName','Analytic')
    end
    plot(times,ave_x_store(:,model_setup),'LineWidth',1.5,'DisplayName',models(model_setup),'Color',setup_colours(model_setup,:))
    xlabel('Time','FontSize',12)
    ylabel('x','FontSize',12)
    legend('FontSize',12)

end

figure(5)
scatter(log(deg_empirical),deviation_store(:,1))
xlabel('log(Degree)')
ylabel('Deviation')

ave_x_store(:,num_models+1) = x_analytic;
expected_adjacency_matrix(:,:,num_models+1) = expected_adjacency_matrix(:,:,3);


%%%%%%%%%% Dynamical simialrity (using x)


% Create dynamical similarity graph
for setup_i=1:num_models+1
    for setup_j=setup_i+1:num_models+1
        dynamical_distance_x(setup_i,setup_j) = norm(ave_x_store(:,setup_i)-ave_x_store(:,setup_j));
        dynamical_distance_x(setup_j,setup_i) = dynamical_distance_x(setup_i,setup_j);
    end
end
dynamicalSimilarity_x = 1./(1 + dynamical_distance_x);
dynamicalSimilarity_x_Graph = graph(dynamicalSimilarity_x - eye(num_models+1));

% Plot dynamical similarity graph
figure (3)
subplot(1,2,1)
pl1 = plot(dynamicalSimilarity_x_Graph,'Layout','circle');
pl1.NodeLabel = {};
pl1.LineWidth = 3.5.*dynamicalSimilarity_x_Graph.Edges.Weight/max(dynamicalSimilarity_x_Graph.Edges.Weight);
pl1.NodeColor = 'r';
pl1.MarkerSize = 10;
pl1.EdgeColor = [0.8 0.8 0.8].*(1 - dynamicalSimilarity_x_Graph.Edges.Weight/max(dynamicalSimilarity_x_Graph.Edges.Weight));
x_adjust = [-.22,.03,.03,-.28,-0.3,-0.2, 0];
y_adjust = [.1,.03,.03,.10,-0.22,-0.19, -0.2];
text(pl1.XData+x_adjust, pl1.YData+y_adjust ,models, ...
    'VerticalAlignment','Bottom',...
    'HorizontalAlignment', 'left',...
    'FontSize', 12)
title(sprintf('Dynamical Similarity Network (t = %.0f)',maxT))


%%%%%%%%%% Network simialrity


% Create network similarity graph
for setup_i=1:num_models
    for setup_j=setup_i+1:num_models
        network_distance(setup_i,setup_j) = calculate_DeltaCon_distance(expected_adjacency_matrix(:,:,setup_i),expected_adjacency_matrix(:,:,setup_j));
        network_distance(setup_j,setup_i) = network_distance(setup_i,setup_j);
    end
end
networkSimilarity = 1./(1 + network_distance);
networkSimilarityGraph = graph(networkSimilarity - eye(num_models));

% Plot network similarity graph
figure (3)
subplot(1,2,2)
pl2 = plot(networkSimilarityGraph,'Layout','circle');
pl2.NodeLabel = {};
pl2.LineWidth = 0.5 + 3.5.*networkSimilarityGraph.Edges.Weight/max(networkSimilarityGraph.Edges.Weight);
pl2.NodeColor = 'r';
pl2.MarkerSize = 10;
pl2.EdgeColor = [0.8 0.8 0.8].*(1 - networkSimilarityGraph.Edges.Weight/max(networkSimilarityGraph.Edges.Weight));
x_adjust = [-.13,.03,-.4,-.28,-0.3,-0.2];
y_adjust = [.05,.03,.05,.06,-0.18,-0.18];
text(pl2.XData+x_adjust, pl2.YData+y_adjust ,models2, ...
    'VerticalAlignment','Bottom',...
    'HorizontalAlignment', 'left',...
    'FontSize', 12)
title('Network Similarity Network')


%%%%%%%%%% Plot expected adjacency matrices


figure(4)
tl5 = tiledlayout(1,num_models,"TileIndexing","columnmajor","TileSpacing","compact");
top = max(max(max(expected_adjacency_matrix)));
bot = min(min(min(expected_adjacency_matrix)));
for model_setup=1:num_models
    nexttile
    imagesc(expected_adjacency_matrix(:,:,model_setup))
    caxis manual
    caxis([bot top]);
    title(models(model_setup),'FontSize',12)
end
colormap lines