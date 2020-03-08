function original_obstacle_layers = GenerateOriginalObstacleLayers()
global xyt_graph_search_ obstacles_
original_obstacle_layers = cell(1, xyt_graph_search_.num_nodes_t);
ang = linspace(0, 2 * pi, 200);
cos_ang = cos(ang);
sin_ang = sin(ang);
for ii = 1 : xyt_graph_search_.num_nodes_t
    x_obs = [];
    y_obs = [];
    for jj = 1 : size(obstacles_, 2)
        obs = obstacles_{1, jj};
        xc = obs{1, ii}.x; yc = obs{1, ii}.y; R = obs{1, ii}.radius;
        x = xc + R .* cos_ang;
        y = yc + R .* sin_ang;
        x_obs = [x_obs, x];
        y_obs = [y_obs, y];
    end
    temp.x = x_obs; temp.y = y_obs;
    original_obstacle_layers{1,ii} = temp;
end
end