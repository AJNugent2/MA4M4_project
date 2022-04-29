% Gillespie simulation for network model

function [G,states] = simulate_exclusion_process(C,community,G,timepoints)

N = sum(C);

% Store initial state
S = community';

% Get edges from network
[sOut,tOut] = findedge(G);

% Set up time
maxT = max(timepoints);
t=0;
idx = 1;
nextTimepoint = timepoints(idx);

% Set up empty array for outputs
states = zeros(length(timepoints),N);

while t<maxT

    % Calculate the rate of each edges
    % The rate along an edge is 1 if the states are different and 0 if they
    % are the same
    rates = (S(sOut)~=S(tOut))'.*G.Edges.Weight;
    rates = rates./N;
    totalRate = sum(rates);

    if totalRate>0

        r1 = rand;
        dt = -log(r1)/totalRate;
        t = t +dt;

        r2 = rand*totalRate;
        edgeChosen = find(cumsum(rates)>=r2,1);
        node1 = sOut(edgeChosen);
        node2 = tOut(edgeChosen);
        S_old = S;
        S(node1) = S_old(node2);
        S(node2) = S_old(node1);

    else

        t = maxT;
        S_old = S;

    end

    % Decide if we need to store this. We need to store if we have gone
    % beyond the next timepoint. We may have crossed multiple timepoints in
    % a single jump.
    while t>nextTimepoint
        states(idx,:) = S_old;
        idx = idx+1;
        if idx<=length(timepoints)
            nextTimepoint = timepoints(idx);
        else
            nextTimepoint = Inf;
        end
    end
    
end

end