function [G,A] = create_ER_network(N,p)

    A = zeros(N,N);

    for i=1:N
        for j=i+1:N
            A(i,j) = (rand<p);
            A(j,i) = A(i,j);
        end
    end

    G = graph(A);

end