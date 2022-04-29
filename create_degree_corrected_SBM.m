function [G,A] = create_degree_corrected_SBM(N,expectedAdjacency)

    A = zeros(N,N);

    for i=1:N
        for j=i+1:N

            A(i,j) = poissrnd(expectedAdjacency(i,j));
            A(j,i) = A(i,j);

        end
    end

    G = graph(A);

end