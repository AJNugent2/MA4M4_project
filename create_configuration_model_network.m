
% Creating a network from deg 

function [G,A] = create_configuration_model_network(deg)

    G = graph();
    G = addnode(G,length(deg));

    % Create a vector of unassigned edges
    X = deg;
    % We need to have two nodes available to connect
    while sum(X>0)>=2
        % Choose a node, with the choice weighted by number of unassigned
        % edges
        r1 = rand*sum(X);
        endpoint1 = find(cumsum(X)>r1,1);
        % Next need to choose the second endpoint of the edge. This must be
        % different from the first endpoint, so temporarily remove that
        % node.
        X_temp = X;
        X_temp(endpoint1) = 0;
        % Now choose the second endpoint from this adapted list 
        r2 = rand*sum(X_temp);
        endpoint2 = find(cumsum(X_temp)>r2,1);
        % Now create an edge
        G = addedge(G,endpoint1,endpoint2,1);
        % Now remove the number of unassigned edges at the endpoints by 1
        X(endpoint1) = X(endpoint1)-1;
        X(endpoint2) = X(endpoint2)-1;
    end

    A = adjacency(G);

end