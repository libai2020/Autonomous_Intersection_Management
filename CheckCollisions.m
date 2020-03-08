function is_collided = CheckCollisions()
is_collided = 1;
global NLP_
Nfe = NLP_.Nfe;
load x.txt
load y.txt
load theta.txt
Nv = length(x) / Nfe;
x = reshape(x',Nfe,Nv)';
y = reshape(y',Nfe,Nv)';
theta = reshape(theta',Nfe,Nv)';
global NLP_
number_of_frame = NLP_.Nfe;
Nv = size(x,1);
global vehicle_geometrics_
lw = vehicle_geometrics_.wheelbase;
hlb = vehicle_geometrics_.width / 2;
global obstacles_ xyt_graph_search_
Nobs = size(obstacles_,2);
index_list = round(linspace(1, xyt_graph_search_.num_nodes_t, number_of_frame));
for ii = 1 : number_of_frame
    for kk = 1 : Nobs
        cur_obs = obstacles_{1,kk};
        xc = cur_obs{1,index_list(ii)}.x;
        yc = cur_obs{1,index_list(ii)}.y;
        rc = cur_obs{1,index_list(ii)}.radius;
        
        for jj = 1 : Nv
            xv = x(jj,ii);
            yv = y(jj,ii);
            tv = theta(jj,ii);
            
            xr = xv + vehicle_geometrics_.r2p * cos(tv);
            yr = yv + vehicle_geometrics_.r2p * sin(tv);
            xf = xv + vehicle_geometrics_.f2p * cos(tv);
            yf = yv + vehicle_geometrics_.f2p * sin(tv);
            
            if (norm([xr - xc, yr - yc]) - rc - vehicle_geometrics_.radius < -0.001)
                return;
            end
            if (norm([xf - xc, yf - yc]) - rc - vehicle_geometrics_.radius < -0.001)
                return;
            end
        end
    end
    
    for jj = 1 : (Nv - 1)
        for kk = (jj + 1) : Nv
            xv1 = x(jj,ii);
            yv1 = y(jj,ii);
            tv1 = theta(jj,ii);
            xv2 = x(kk,ii);
            yv2 = y(kk,ii);
            tv2 = theta(kk,ii);
            if (IsConfig1CollidingWithConfig2([xv1, yv1, tv1], [xv2, yv2, tv2]))
                return;
            end
        end
    end
end
is_collided = 0;
end

function is_colliding = IsConfig1CollidingWithConfig2(vec1, vec2)
global vehicle_geometrics_
xr = [vec1(1), vec2(1)] +  vehicle_geometrics_.r2p * cos([vec1(3), vec2(3)]);
yr = [vec1(2), vec2(2)] +  vehicle_geometrics_.r2p * sin([vec1(3), vec2(3)]);
xf = [vec1(1), vec2(1)] +  vehicle_geometrics_.f2p * cos([vec1(3), vec2(3)]);
yf = [vec1(2), vec2(2)] +  vehicle_geometrics_.f2p * sin([vec1(3), vec2(3)]);
P11 = [xr(1), yr(1)];
P12 = [xf(1), yf(1)];
P21 = [xr(2), yr(2)];
P22 = [xf(2), yf(2)];
distance = min([norm(P11 - P21), norm(P11 - P22), norm(P12 - P21), norm(P12 - P22)]) - 2 * vehicle_geometrics_.radius;

if (distance <= -0.001)
    is_colliding = 1;
else
    is_colliding = 0;
end
end