
% Creating a network from degWithin and degBetween

function [G,A] = create_community_configuration_network(degWithin,degBetween,C1)

    G = graph();
    G = addnode(G,length(degWithin));

    % Create edges within C1
    % Create a vector of unassigned edges
    X1 = degWithin(1:C1);
    % We need to have two nodes available to connect
    while sum(X1>0)>=2
        % Choose a node, with the choice weighted by number of unassigned
        % edges
        r1 = rand*sum(X1);
        endpoint1 = find(cumsum(X1)>r1,1);
        % Next need to choose the second endpoint of the edge. This must be
        % different from the first endpoint, so temporarily remove that
        % node.
        X1_temp = X1;
        X1_temp(endpoint1) = 0;
        % Now choose the second endpoint from this adapted list 
        r2 = rand*sum(X1_temp);
        endpoint2 = find(cumsum(X1_temp)>r2,1);
        % Now create an edge
        G = addedge(G,endpoint1,endpoint2,1);
        % Now remove the number of unassigned edges at the endpoints by 1
        X1(endpoint1) = X1(endpoint1)-1;
        X1(endpoint2) = X1(endpoint2)-1;
    end

    % Create edges within C2
    % Create a vector of unassigned edges
    X2 = degWithin(1+C1:end);
    % We need to have two nodes available to connect
    while sum(X2>0)>=2
        % Choose a node, with the choice weighted by number of unassigned
        % edges
        r1 = rand*sum(X2);
        endpoint1 = find(cumsum(X2)>r1,1);
        % Next need to choose the second endpoint of the edge. This must be
        % different from the first endpoint, so temporarily remove that
        % node.
        X2_temp = X2;
        X2_temp(endpoint1) = 0;
        % Now choose the second endpoint from this adapted list 
        r2 = rand*sum(X2_temp);
        endpoint2 = find(cumsum(X2_temp)>r2,1);
        % Now create an edge
        G = addedge(G,C1+endpoint1,C1+endpoint2,1);
        % Now remove the number of unassigned edges at the endpoints by 1
        X2(endpoint1) = X2(endpoint1)-1;
        X2(endpoint2) = X2(endpoint2)-1;
    end

    % Create edges between C1 and C2
    % Create a vector of unassigned edges
    X1to2 = degBetween(1:C1);
    X2to1 = degBetween(1+C1:end);
    % We need to have at least nodes available to connect in each community
    while sum(X1to2)>0 && sum(X2to1)>0
        % Choose a node in C1, with the choice weighted by number of unassigned
        % edges
        r1 = rand*sum(X1to2);
        C1endpoint = find(cumsum(X1to2)>r1,1);
        % Choose a node in C2 to connect to
        r2 = rand*sum(X2to1);
        C2endpoint = find(cumsum(X2to1)>r2,1);
        % Now create an edge
        G = addedge(G,C1endpoint,C1+C2endpoint,1);
        % Now remove the number of unassigned edges at the endpoints by 1
        X1to2(C1endpoint) = X1to2(C1endpoint)-1;
        X2to1(C2endpoint) = X2to1(C2endpoint)-1;
    end

    A = adjacency(G);

end