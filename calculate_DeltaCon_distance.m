function d = calculate_DeltaCon_distance(A1,A2)

    N = length(A1);

    d1 = sum(A1);
    D1 = eye(N).*d1;

    d2 = sum(A2);
    D2 = eye(N).*d2;

    e1 = 1/(1 + max(d1));
    e2 = 1/(1 + max(d2));

    S1 = inv(eye(N) + D1.*e1^2 - A1.*e1);
    S2 = inv(eye(N) + D2.*e2^2 - A2.*e2);

    toSum = zeros(N,N);
    for i=1:N
        for j=1:N
            toSum(i,j) = (sqrt(S1(i,j)) - sqrt(S2(i,j)))^2;
        end
    end

    d = sqrt(sum(sum(toSum)));

end