clear G_empirical community

% Load data on polblogs network from gml file
S = fileread('polblogs.gml');
% Create graph from gml data
%nodes = regexp(S, 'node.*?id (?<id>\d+).*?label\s*"(?<label>[^"]*)"', 'names');
nodes = regexp(S, 'node.*?id (?<id>\d+).*?label\s*"(?<label>[^"]*)".*?value (?<value>\d+)', 'names');
edges = regexp(S, 'edge.*?source\s*(?<source>\d+).*?target\s*(?<target>\d+)', 'names');
all_ids = {nodes.id};
all_names = {nodes.label};
all_values = {nodes.value};
all_sources = {edges.source};
all_targets = {edges.target};
[source_found, s] = ismember(all_sources, all_ids);
[target_found, t] = ismember(all_targets, all_ids);
EdgeTable = table([s.', t.'], ones(length(s),1), 'VariableNames', {'EndNodes' 'Weight'});
%NodeTable = table(all_names.', 'VariableNames',{'Name'});
NodeTable = table(all_values.', 'VariableNames',{'value'});
G_empirical = graph(EdgeTable,NodeTable);

% Restrict to largest connected component of the network 
[bin,binsize] = conncomp(G_empirical);
idxs = binsize(bin)==max(binsize);
G_empirical = subgraph(G_empirical,idxs);

% Remove any self- or multi- edges
G_empirical = simplify(G_empirical);
G_empirical.Nodes.Idx = (1:1:height(G_empirical.Nodes))';

C0 = 586;
N = 1222;

community = zeros(height(G_empirical.Nodes),1);
community(1:C0) = 0;
community(1+C0:N) = 1;

% Add this property to the nodes of G 
G_empirical.Nodes.Comm = community;


% Sort nodes and community vector in community order
[~,order] = sort(community,'ascend');
G_empirical = reordernodes(G_empirical,order);
community = community(order);

% Some useful network features
C0 = sum(community==0);
C1 = sum(community==1);
C = [C0,C1];
numC = length(C);
N = sum(C);
c1=C0/N;
c2=C1/N;

% Sort nodes within communities by degree order
deg_empirical = degree(G_empirical);
[~,orderC0] = sort(deg_empirical(1:C0),'descend');
[~,orderC1] = sort(deg_empirical(1+C0:end),'descend');
orderC1 = orderC1 + C0;
G_empirical = reordernodes(G_empirical,[orderC0' orderC1']);
deg_empirical = degree(G_empirical);
G_empirical.Nodes.Degree = deg_empirical;

A_empirical = adjacency(G_empirical);

% Calculate probability of choosing a neighbour in the other community
[m,deg_w_empirical,deg_b_empirical] = get_network_mixing_rate(G_empirical,C0,C1);

colours = [0.00 0.45 0.74;0.8500 0.3250 0.098;0.494 0.184 0.5560;0.4660 0.6740 0.1880;0 0 0;0.1 0.4 0.3 ];


% The following is used to setup the degree-corrected SBM: 

m_01 = sum(deg_b_empirical)/2;
m_00 = sum(deg_w_empirical(1:C0));
m_11 = sum(deg_w_empirical(1+C0:end));

theta = zeros(1,N);
theta(1:C0) = deg_empirical(1:C0)/sum(deg_empirical(1:C0));
theta(1+C0:end) = deg_empirical(1+C0:end)/sum(deg_empirical(1+C0:end));

deg_corr_SBM_expected_adjacency = zeros(N,N);
for i=1:C0
    for j=i+1:C0
        deg_corr_SBM_expected_adjacency(i,j) = theta(i)*theta(j)*m_00;
        deg_corr_SBM_expected_adjacency(j,i) = deg_corr_SBM_expected_adjacency(i,j);
    end
    for j=C0+1:N
        deg_corr_SBM_expected_adjacency(i,j) = theta(i)*theta(j)*m_01;
        deg_corr_SBM_expected_adjacency(j,i) = deg_corr_SBM_expected_adjacency(i,j);
    end
end
for i=1+C0:N
    for j=i+1:N
        deg_corr_SBM_expected_adjacency(i,j) = theta(i)*theta(j)*m_11;
        deg_corr_SBM_expected_adjacency(j,i) = deg_corr_SBM_expected_adjacency(i,j);
    end
end