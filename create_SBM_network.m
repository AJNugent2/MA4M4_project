function [G,A] = create_SBM_network(C,p_in,p_out)

    numC = length(C); % Number of communities
    N = sum(C);       % Total number of nodes

    A = zeros(N,N);

    for c=1:numC

        for i=1+sum(C(1:c-1)):sum(C(1:c))
            % Nodes within communities
            for j=i+1:sum(C(1:c))
                A(i,j) = (rand<p_in(c));
                A(j,i) = A(i,j);
            end
            % Nodes between communities
            for j=sum(C(1:c)):N
                A(i,j) = (rand<p_out);
                A(j,i) = A(i,j);
            end
        end

    end

    G = graph(A);

end