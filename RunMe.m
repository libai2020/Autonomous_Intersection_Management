%  MATLAB Source Codes for the book "Cooperative Decision and Planning for
%  Connected and Automated Vehicles" published by Mechanical Industry Press
%  in 2020.
% 《智能网联汽车协同决策与规划技术》书籍配套代码
%  Copyright (C) 2020 Bai Li
%  2020.03.09
% ==============================================================================
%  第五章 基于X-Y-T A*算法热启动以及增量初始化辅助求解的无信号灯路口多车协同轨迹规划方法
% ==============================================================================
%  备注1. 如果内存足够大，可以将Nv设置得更大；
%  备注2. 由于MA57不稳定，偶尔出现线性求解器崩溃情况，此问题未来亟待修复，可在ipopt.opt文件中将ma57禁用；
%  备注3. 使用该部分代码二次开发的读者在著作中宜引用以下直接相关的文献：
%  a) Bai Li, Youmin Zhang, Ning Jia, and Xiaoyan Peng, "Autonomous
%  Intersection Management over Continuous Space: A Microscopic and Precise
%  Solution by Computational Optimal Control", IFAC 2020, accepted, 2020.
%  b) Bai Li, and Youmin Zhang, "Fault-tolerant cooperative motion planning
%  of connected and automated vehicles at a signal-free and lane-free
%  intersection." IFAC-PapersOnLine 51.24 (2018): 60-67.
%  c) Li, B., Zhang, Y. M., Zhang, Y., Jia, N., & Ge, Y. (2018).
%  Near-optimal online motion planning of connected and automated vehicles
%  at a signal-free and lane-free intersection. In 2018 IEEE Intelligent
%  Vehicles Symposium (IV) (pp. 1432-1437). IEEE.
%  d) Li, B., Jia, N., Li, P., Ran, X., & Li, Y. (2019). Incrementally
%  constrained dynamic optimization: A computational framework for lane
%  change motion planning of connected and automated vehicles. Journal of
%  Intelligent Transportation Systems, 23(6), 557-568.
% ==============================================================================
clear; close all; clc;
% % Parameters
global vehicle_geometrics_
vehicle_geometrics_.wheelbase = 2.8;
vehicle_geometrics_.front_hang = 0.96;
vehicle_geometrics_.rear_hang = 0.929;
vehicle_geometrics_.width = 1.942;
vehicle_geometrics_.length = vehicle_geometrics_.wheelbase + vehicle_geometrics_.front_hang + vehicle_geometrics_.rear_hang;
vehicle_geometrics_.radius = hypot(0.25 * vehicle_geometrics_.length, 0.5 * vehicle_geometrics_.width);
vehicle_geometrics_.r2p = 0.25 * vehicle_geometrics_.length - vehicle_geometrics_.rear_hang;
vehicle_geometrics_.f2p = 0.75 * vehicle_geometrics_.length - vehicle_geometrics_.rear_hang;
global vehicle_kinematics_
vehicle_kinematics_.v_max = 25;
vehicle_kinematics_.v_common = 20;
vehicle_kinematics_.a_max = 2.0;
vehicle_kinematics_.phy_max = 0.7;
vehicle_kinematics_.w_max = 0.3;
global environment_scale_
environment_scale_.x_min = -250;
environment_scale_.x_max = 250;
environment_scale_.y_min = -250;
environment_scale_.y_max = 250;
environment_scale_.x_scale = environment_scale_.x_max - environment_scale_.x_min;
environment_scale_.y_scale = environment_scale_.y_max - environment_scale_.y_min;
global intersection
intersection.road_width = 12;
intersection.num_lanes = 3;
global xyt_graph_search_
xyt_graph_search_.t_max = 10;
xyt_graph_search_.num_nodes_t = 200;
xyt_graph_search_.resolution_t = xyt_graph_search_.t_max / xyt_graph_search_.num_nodes_t;
xyt_graph_search_.resolution_x = 1; xyt_graph_search_.resolution_y = 1;
xyt_graph_search_.num_nodes_x = ceil(environment_scale_.x_scale / xyt_graph_search_.resolution_x) + 1;
xyt_graph_search_.num_nodes_y = ceil(environment_scale_.y_scale / xyt_graph_search_.resolution_y) + 1;
xyt_graph_search_.multiplier_for_H = 2.0;
xyt_graph_search_.weight_for_time = 5;
xyt_graph_search_.max_cycle = 2000;
xyt_graph_search_.pelnaty_for_halt = 20;
ang = linspace(0, 2 * pi, 200);
xyt_graph_search_.Rcos_ang = vehicle_geometrics_.radius .* cos(ang);
xyt_graph_search_.Rsin_ang = vehicle_geometrics_.radius .* sin(ang);
global NLP_
NLP_.Nfe = 100;
NLP_.tf = xyt_graph_search_.t_max;
global boundary_configs_ obstacles_
Nv = 24;
[boundary_configs_, obstacles_] = GenerateAIM(Nv);
[index_flag, ranklist] = SpecifyRanklist(Nv);

global original_obstacle_layers_
original_obstacle_layers_ = GenerateOriginalObstacleLayers();
backup_original_obstacle_layers = original_obstacle_layers_;

DrawIntersection();
decision_x = zeros(Nv, NLP_.Nfe);
decision_y = zeros(Nv, NLP_.Nfe);
decision_theta = zeros(Nv, NLP_.Nfe);
for ii = 1 : Nv
    iv = ranklist(ii);
    [x, y, theta] = SearchTrajectoryInXYTGraph(boundary_configs_{1,iv});
    plot(x,y); hold on; drawnow;
    decision_x(iv, 1 : NLP_.Nfe) = ResampleProfile(x, NLP_.Nfe);
    decision_y(iv, 1 : NLP_.Nfe) = ResampleProfile(y, NLP_.Nfe);
    decision_theta(iv, 1 : NLP_.Nfe) = ResampleProfile(theta, NLP_.Nfe);
end
% ShowDynamicsInitialGuess();

WriteBoundaryValuesAndBasicParams(Nv, NLP_.Nfe);
SpecifyInitialGuess(decision_x, decision_y, decision_theta);
!ampl rr0.run
load opti_flag.txt
if (~opti_flag)
    error 'Solution_Failed.'
end

for ii = 1 : 2 : NLP_.Nfe
    delete('TimeLine');
    fid = fopen('TimeLine', 'w');
    fprintf(fid, '1 %g\r\n', ii);
    fclose(fid);
    !ampl rr.run
    load opti_flag.txt
    if (opti_flag && (~CheckCollisions()))
        break;
    end
end
if (opti_flag && (~CheckCollisions()))
    disp('Solution_Succeeded.');
    ShowDynamicResults();
else
    error 'Solution_Failed.'
end