function Rc = coarsen(R, m)
    mc = (m - 1) / 2;
    R = reshape(R, m, m);  % Convert to 2D

    Rc = zeros(mc, mc);

    for i = 1:mc
        for j = 1:mc
            ii = 2 * i;
            jj = 2 * j;

            % Ensure indices are within bounds
            val = 0;

            for di = -1:1
                for dj = -1:1
                    weight = 1;
                    if abs(di) + abs(dj) == 1
                        weight = 2;
                    elseif di == 0 && dj == 0
                        weight = 4;
                    end

                    ni = ii + di;
                    nj = jj + dj;

                    if ni >= 1 && ni <= m && nj >= 1 && nj <= m
                        val = val + weight * R(ni, nj);
                    end
                end
            end

            Rc(i, j) = val / 16;
        end
    end

    Rc = Rc(:);  % Flatten
end
