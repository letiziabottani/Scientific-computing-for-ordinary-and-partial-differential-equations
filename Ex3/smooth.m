function Unew = smooth(U, omega, m, F)
    h = 1/(m+1);
    U = reshape(U, m, m);  % reshape into 2D grid
    F = reshape(F, m, m);
    Unew = U;

    for i = 2:m-1
        for j = 2:m-1
            Unew(i,j) = (1 - omega) * U(i,j) + omega * 0.25 * ...
                (U(i+1,j) + U(i-1,j) + U(i,j+1) + U(i,j-1) - h^2 * F(i,j));
        end
    end

    Unew = reshape(Unew, m*m, 1);  % flatten back to vector
end
