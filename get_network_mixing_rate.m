
function [m,deg_within,deg_between] = get_network_mixing_rate(G,C0,C1)

    N = C0+C1;
    
    % Create empty arrays for within- and between-degree
    deg_within = zeros(1,N);
    deg_between = zeros(1,N);
    
    % Loop over nodes in C0
    for i=1:C0
        deg_within(i) = full(sum(edgecount(G,i,1:C0)));
        deg_between(i) = full(sum(edgecount(G,i,C0+1:N)));
    end

    % Loop over nodes in C1
    for i=1:C1
        j = i+C0;
        deg_within(j) = full(sum(edgecount(G,j,C0+1:N)));
        deg_between(j) = full(sum(edgecount(G,j,1:C0)));
    end

    % Calculate the number of edges between communities
    Mb = sum(deg_between)/2;
    % Calculate m 
    m = (N/C0)*(N/C1)*(Mb)/N;

end