function R = interpolate(Rc, m)
    mc = (m - 1) / 2;
    Rc = reshape(Rc, mc, mc);
    R = zeros(m, m);

    for i = 1:mc
        for j = 1:mc
            ii = 2 * i;
            jj = 2 * j;
            value = Rc(i, j);

            for di = -1:1
                for dj = -1:1
                    ni = ii + di;
                    nj = jj + dj;

                    if ni >= 1 && ni <= m && nj >= 1 && nj <= m
                        weight = 1;
                        if abs(di) + abs(dj) == 1
                            weight = 2;
                        elseif di == 0 && dj == 0
                            weight = 4;
                        end

                        R(ni, nj) = R(ni, nj) + (weight / 4) * value;
                    end
                end
            end
        end
    end

    R = R(:);
end
