function [index_flag, ranklist] = SpecifyRanklist(Nv)
distance_to_origin = zeros(1, Nv);
index_flag = zeros(1, Nv);
global A1 A2 A3 B1 B2 B3 C1 C2 C3 D1 D2 D3 boundary_configs_
for ii = 1 : (A1 + A2 + A3)
    distance_to_origin(ii) = 10 * abs(boundary_configs_{1,ii}.x0) + abs(boundary_configs_{1,ii}.y0);
    index_flag(ii) = 1;
end

for ii = (A1 + A2 + A3 + 1) : (A1 + A2 + A3 + B1 + B2 + B3)
    distance_to_origin(ii) = abs(boundary_configs_{1,ii}.x0) + 10 * abs(boundary_configs_{1,ii}.y0);
    index_flag(ii) = 2;
end

for ii = (A1 + A2 + A3 + B1 + B2 + B3 + 1) : (A1 + A2 + A3 + B1 + B2 + B3 + C1 + C2 + C3)
    distance_to_origin(ii) = 10 * abs(boundary_configs_{1,ii}.x0) + abs(boundary_configs_{1,ii}.y0);
    index_flag(ii) = 3;
end

for ii = (A1 + A2 + A3 + B1 + B2 + B3 + C1 + C2 + C3 + 1) : Nv
    distance_to_origin(ii) = abs(boundary_configs_{1,ii}.x0) + 10 * abs(boundary_configs_{1,ii}.y0);
    index_flag(ii) = 4;
end

[~, ranklist] = sort(distance_to_origin,'descend');
end