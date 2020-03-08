param cav_index_start{i in {1..12}};
param cav_index_end{i in {1..12}};
set A1 := {cav_index_start[1]..cav_index_end[1]};
set A2 := {cav_index_start[2]..cav_index_end[2]};
set A3 := {cav_index_start[3]..cav_index_end[3]};
set B1 := {cav_index_start[4]..cav_index_end[4]};
set B2 := {cav_index_start[5]..cav_index_end[5]};
set B3 := {cav_index_start[6]..cav_index_end[6]};
set C1 := {cav_index_start[7]..cav_index_end[7]};
set C2 := {cav_index_start[8]..cav_index_end[8]};
set C3 := {cav_index_start[9]..cav_index_end[9]};
set D1 := {cav_index_start[10]..cav_index_end[10]};
set D2 := {cav_index_start[11]..cav_index_end[11]};
set D3 := {cav_index_start[12]..cav_index_end[12]};

param BP{i in {1..12}};
param wheelbase == BP[1];
param FourSqR == 4 * ((BP[2])^2);
param R == BP[2];
param r2p == BP[3];
param f2p == BP[4];
param vmax == BP[5];
param amax == BP[6];
param phymax == BP[7];
param wmax == BP[8];
param tf == BP[9];
param Nv == BP[10];
param Nfe == BP[11];
param vcommon == BP[12];
param hi = tf / Nfe;

param BV{i in {1..Nv}, j in {1..6}};
param TimeLine{i in {1..1}};
param Nt == TimeLine[1];
set I := {1..Nt};
set IPmnt := {1..Nfe};
set IP := {1..Nfe};

var x{i in {1..Nv}, j in {1..Nfe}};
var y{i in {1..Nv}, j in {1..Nfe}};
var xf{i in {1..Nv}, j in {1..Nfe}};
var yf{i in {1..Nv}, j in {1..Nfe}};
var xr{i in {1..Nv}, j in {1..Nfe}};
var yr{i in {1..Nv}, j in {1..Nfe}};
var theta{i in {1..Nv}, j in {1..Nfe}};
var v{i in {1..Nv}, j in {1..Nfe}};
var a{i in {1..Nv}, j in {1..Nfe}};
var phy{i in {1..Nv}, j in {1..Nfe}};
var w{i in {1..Nv}, j in {1..Nfe}};

minimize cost_function:
hi * sum{i in {1..Nv}, j in I}(a[i,j]^2 + w[i,j]^2) + sum{k in A1}(x[k,Nfe]) + sum{k in A2}(-y[k,Nfe]) + sum{k in A3}(-x[k,Nfe]) + sum{k in B1}(y[k,Nfe]) + sum{k in B2}(x[k,Nfe]) + sum{k in B3}(-y[k,Nfe]) + sum{k in C1}(-x[k,Nfe]) + sum{k in C2}(y[k,Nfe]) + sum{k in C3}(x[k,Nfe]) + sum{k in D1}(-y[k,Nfe]) + sum{k in D2}(-x[k,Nfe]) + sum{k in D3}(y[k,Nfe]);


s.t. DIFF_dxdt {k in {1..Nv}, i in {2..Nfe}}:
x[k,i] = x[k,i-1] + hi * v[k,i] * cos(theta[k,i]);
s.t. DIFF_dydt {k in {1..Nv}, i in {2..Nfe}}:
y[k,i] = y[k,i-1] + hi * v[k,i] * sin(theta[k,i]);
s.t. DIFF_dvdt {k in {1..Nv}, i in {2..Nfe}}:
v[k,i] = v[k,i-1] + hi * a[k,i];
s.t. DIFF_dthetadt {k in {1..Nv}, i in {2..Nfe}}:
theta[k,i] = theta[k,i-1] + hi * tan(phy[k,i]) * v[k,i] / wheelbase;
s.t. DIFF_dphydt {k in {1..Nv}, i in {2..Nfe}}:
phy[k,i] = phy[k,i-1] + hi * w[k,i];

s.t. RELATIONSHIP_XF {k in {1..Nv}, i in IPmnt}:
xf[k,i] = x[k,i] + f2p * cos(theta[k,i]);
s.t. RELATIONSHIP_YF {k in {1..Nv}, i in IPmnt}:
yf[k,i] = y[k,i] + f2p * sin(theta[k,i]);
s.t. RELATIONSHIP_XR {k in {1..Nv}, i in IPmnt}:
xr[k,i] = x[k,i] + r2p * cos(theta[k,i]);
s.t. RELATIONSHIP_YR {k in {1..Nv}, i in IPmnt}:
yr[k,i] = y[k,i] + r2p * sin(theta[k,i]);

s.t. Bonds_phy {k in {1..Nv}, i in {1..Nfe}}:
-phymax <= phy[k,i] <= phymax;
s.t. Bonds_a {k in {1..Nv}, i in {1..Nfe}}:
-amax <= a[k,i] <= amax;
s.t. Bonds_v {k in {1..Nv}, i in {1..Nfe}}:
0 <= v[k,i] <= vmax;
s.t. Bonds_w {k in {1..Nv}, i in {1..Nfe}}:
-wmax <= w[k,i] <= wmax;

### D1 ###
s.t. L1argeCircle_for_D1 {i in D1, j in IP}:
((xr[i,j]+100)^2 + (yr[i,j]-100)^2) / 8014.224 >= 1;
s.t. M1iddleCircle_for_D1 {i in D1, j in IP}:
((xr[i,j]+56)^2 + (yr[i,j]-56)^2) / 2072.271 >= 1;
s.t. S1mallCircle_for_D1 {i in D1, j in IP}:
((xr[i,j]+34)^2 + (yr[i,j]-34)^2) / 553.294 >= 1;
s.t. L2argeCircle_for_D1 {i in D1, j in IP}:
((xf[i,j]+100)^2 + (yf[i,j]-100)^2) / 8014.224 >= 1;
s.t. M2iddleCircle_for_D1 {i in D1, j in IP}:
((xf[i,j]+56)^2 + (yf[i,j]-56)^2) /  2072.271 >= 1;
s.t. S2mallCircle_for_D1 {i in D1, j in IP}:
((xf[i,j]+34)^2 + (yf[i,j]-34)^2) / 553.294 >= 1;
s.t. limit_xr_D1 {i in D1, j in IP}:
xr[i,j] + R <= 12;
s.t. limit_xf_D1 {i in D1, j in IP}:
xf[i,j] + R <= 12;
s.t. limit_yr_D1 {i in D1, j in IP}:
yr[i,j] - R >= -12;
s.t. limit_yf_D1 {i in D1, j in IP}:
yf[i,j] - R >= -12;
s.t. terminal_x_D1 {m in D1}:
yr[m,Nfe] - R >= 12;
s.t. terminal_xr_D1 {m in D1}:
xr[m,Nfe] - R >= 0;
s.t. terminal_xf_D1 {m in D1}:
xf[m,Nfe] - R >= 0;
s.t. interval_theta_D1 {i in D1, j in IPmnt}:
-1 <= theta[i,j] <= 2;
### D2 ###
s.t. limit_yr_D2 {i in D2, j in IP}:
-12 + R <= yr[i,j] <= -R;
s.t. limit_yf_D2 {i in D2, j in IP}:
-12 + R <= yf[i,j] <= -R;
s.t. terminal_x_D2 {m in D2}:
xr[m,Nfe] - R >= 12;
s.t. interval_theta_D2 {i in D2, j in IPmnt}:
-1 <= theta[i,j] <= 1;
### D3 ###
s.t. L1argeCircle_for_D3 {i in D3, j in IP}:
((xr[i,j]+100)^2 + (yr[i,j]+100)^2) / 8014.224 >= 1;
s.t. M1iddleCircle_for_D3 {i in D3, j in IP}:
((xr[i,j]+56)^2 + (yr[i,j]+56)^2) /  2072.271 >= 1;
s.t. S1mallCircle_for_D3 {i in D3, j in IP}:
((xr[i,j]+34)^2 + (yr[i,j]+34)^2) / 553.294 >= 1;
s.t. L2argeCircle_for_D3 {i in D3, j in IP}:
((xf[i,j]+100)^2 + (yf[i,j]+100)^2) / 8014.224 >= 1;
s.t. M2iddleCircle_for_D3 {i in D3, j in IP}:
((xf[i,j]+56)^2 + (yf[i,j]+56)^2) /  2072.271 >= 1;
s.t. S2mallCircle_for_D3 {i in D3, j in IP}:
((xf[i,j]+34)^2 + (yf[i,j]+34)^2) / 553.294 >= 1;
s.t. limit_xr_D3 {i in D3, j in IP}:
xr[i,j] + R <= 0;
s.t. limit_xf_D3 {i in D3, j in IP}:
xf[i,j] + R <= 0;
s.t. limit_yr_D3 {i in D3, j in IP}:
yr[i,j] + R <= 0;
s.t. limit_yf_D3 {i in D3, j in IP}:
yf[i,j] + R <= 0;
s.t. terminal_x_D3 {m in D3}:
yr[m,Nfe] + R <= -12;
s.t. interval_theta_D3 {i in D3, j in IPmnt}:
-2 <= theta[i,j] <= 1;
### A1 ###
s.t. L1argeCircle_for_A1 {i in A1, j in IP}:
((xr[i,j]+100)^2 + (yr[i,j]+100)^2) / 8014.224 >= 1;
s.t. M1iddleCircle_for_A1 {i in A1, j in IP}:
((xr[i,j]+56)^2 + (yr[i,j]+56)^2) /  2072.271 >= 1;
s.t. S1mallCircle_for_A1 {i in A1, j in IP}:
((xr[i,j]+34)^2 + (yr[i,j]+34)^2) / 553.294 >= 1;
s.t. L2argeCircle_for_A1 {i in A1, j in IP}:
((xf[i,j]+100)^2 + (yf[i,j]+100)^2) / 8014.224 >= 1;
s.t. M2iddleCircle_for_A1 {i in A1, j in IP}:
((xf[i,j]+56)^2 + (yf[i,j]+56)^2) /  2072.271 >= 1;
s.t. S2mallCircle_for_A1 {i in A1, j in IP}:
((xf[i,j]+34)^2 + (yf[i,j]+34)^2) / 553.294 >= 1;
s.t. limit_xr_A1 {i in A1, j in IP}:
xr[i,j] + R <= 12;
s.t. limit_xf_A1 {i in A1, j in IP}:
xf[i,j] + R <= 12;
s.t. limit_yr_A1 {i in A1, j in IP}:
yr[i,j] + R <= 12;
s.t. limit_yf_A1 {i in A1, j in IP}:
yf[i,j] + R <= 12;
s.t. terminal_x_A1 {m in A1}:
xr[m,Nfe] + R <= -12;
s.t. terminal_yr_A1 {m in A1}:
yr[m,Nfe] - R >= 0;
s.t. terminal_yf_A1 {m in A1}:
yf[m,Nfe] - R >= 0;
s.t. interval_theta_A1 {i in A1, j in IPmnt}:
1 <= theta[i,j] <= 4;
### A2 ###
s.t. limit_xr_A2 {i in A2, j in IP}:
R <= xr[i,j] <= 12 - R;
s.t. limit_xf_A2 {i in A2, j in IP}:
R <= xf[i,j] <= 12 - R;
s.t. terminal_x_A2 {m in A2}:
yr[m,Nfe] - R >= 12;
s.t. interval_theta_A2 {i in A2, j in IPmnt}:
1 <= theta[i,j] <= 2;
### A3 ###
s.t. L1argeCircle_for_A3 {i in A3, j in IP}:
((xr[i,j]-100)^2 + (yr[i,j]+100)^2) / 8014.224 >= 1;
s.t. M1iddleCircle_for_A3 {i in A3, j in IP}:
((xr[i,j]-56)^2 + (yr[i,j]+56)^2) /  2072.271 >= 1;
s.t. S1mallCircle_for_A3 {i in A3, j in IP}:
((xr[i,j]-34)^2 + (yr[i,j]+34)^2) / 553.294 >= 1;
s.t. L2argeCircle_for_A3 {i in A3, j in IP}:
((xf[i,j]-100)^2 + (yf[i,j]+100)^2) / 8014.224 >= 1;
s.t. M2iddleCircle_for_A3 {i in A3, j in IP}:
((xf[i,j]-56)^2 + (yf[i,j]+56)^2) /  2072.271 >= 1;
s.t. S2mallCircle_for_A3 {i in A3, j in IP}:
((xf[i,j]-34)^2 + (yf[i,j]+34)^2) / 553.294 >= 1;
s.t. limit_xr_A3 {i in A3, j in IP}:
xr[i,j] - R >= 0;
s.t. limit_xf_A3 {i in A3, j in IP}:
xf[i,j] - R >= 0;
s.t. limit_yr_A3 {i in A3, j in IP}:
yr[i,j] + R <= 0;
s.t. limit_yf_A3 {i in A3, j in IP}:
yf[i,j] + R <= 0;
s.t. terminal_x_A3 {m in A3}:
xr[m,Nfe] - R >= 12;
s.t. interval_theta_A3 {i in A3, j in IPmnt}:
-1 <= theta[i,j] <= 2;
### B1 ###
s.t. L1argeCircle_for_B1 {i in B1, j in IP}:
((xr[i,j]-100)^2 + (yr[i,j]+100)^2) / 8014.224 >= 1;
s.t. M1iddleCircle_for_B1 {i in B1, j in IP}:
((xr[i,j]-56)^2 + (yr[i,j]+56)^2) /  2072.271 >= 1;
s.t. S1mallCircle_for_B1 {i in B1, j in IP}:
((xr[i,j]-34)^2 + (yr[i,j]+34)^2) / 553.294 >= 1;
s.t. L2argeCircle_for_B1 {i in B1, j in IP}:
((xf[i,j]-100)^2 + (yf[i,j]+100)^2) / 8014.224 >= 1;
s.t. M2iddleCircle_for_B1 {i in B1, j in IP}:
((xf[i,j]-56)^2 + (yf[i,j]+56)^2) /  2072.271 >= 1;
s.t. S2mallCircle_for_B1 {i in B1, j in IP}:
((xf[i,j]-34)^2 + (yf[i,j]+34)^2) / 553.294 >= 1;
s.t. limit_xr_B1 {i in B1, j in IP}:
xr[i,j] - R >= -12;
s.t. limit_xf_B1 {i in B1, j in IP}:
xf[i,j] - R >= -12;
s.t. limit_yr_B1 {i in B1, j in IP}:
yr[i,j] + R <= 12;
s.t. limit_yf_B1 {i in B1, j in IP}:
yf[i,j] + R <= 12;
s.t. terminal_x_B1 {m in B1}:
yr[m,Nfe] + R <= -12;
s.t. terminal_xr_B1 {m in B1}:
xr[m,Nfe] + R <= 0;
s.t. terminal_xf_B1 {m in B1}:
xf[m,Nfe] + R <= 0;
s.t. interval_theta_B1 {i in B1, j in IPmnt}:
2 <= theta[i,j] <= 6;
### B2 ###
s.t. limit_yr_B2 {i in B2, j in IP}:
R <= yr[i,j] <= 12 - R;
s.t. limit_yf_B2 {i in B2, j in IP}:
R <= yf[i,j] <= 12 - R;
s.t. terminal_x_B2 {m in B2}:
xr[m,Nfe] + R <= -12;
s.t. interval_theta_B2 {i in B2, j in IPmnt}:
2 <= theta[i,j] <= 4;
### B3 ###
s.t. L1argeCircle_for_B3 {i in B3, j in IP}:
((xr[i,j]-100)^2 + (yr[i,j]-100)^2) / 8014.224 >= 1;
s.t. M1iddleCircle_for_B3 {i in B3, j in IP}:
((xr[i,j]-56)^2 + (yr[i,j]-56)^2) /  2072.271 >= 1;
s.t. S1mallCircle_for_B3 {i in B3, j in IP}:
((xr[i,j]-34)^2 + (yr[i,j]-34)^2) / 553.294 >= 1;
s.t. L2argeCircle_for_B3 {i in B3, j in IP}:
((xf[i,j]-100)^2 + (yf[i,j]-100)^2) / 8014.224 >= 1;
s.t. M2iddleCircle_for_B3 {i in B3, j in IP}:
((xf[i,j]-56)^2 + (yf[i,j]-56)^2) /  2072.271 >= 1;
s.t. S2mallCircle_for_B3 {i in B3, j in IP}:
((xf[i,j]-34)^2 + (yf[i,j]-34)^2) / 553.294 >= 1;
s.t. limit_xr_B3 {i in B3, j in IP}:
xr[i,j] - R >= 0;
s.t. limit_xf_B3 {i in B3, j in IP}:
xf[i,j] - R >= 0;
s.t. limit_yr_B3 {i in B3, j in IP}:
yr[i,j] - R >= 0;
s.t. limit_yf_B3 {i in B3, j in IP}:
yf[i,j] - R >= 0;
s.t. terminal_x_B3 {m in B3}:
yr[m,Nfe] - R >= 12;
s.t. interval_theta_B3 {i in B3, j in IPmnt}:
1 <= theta[i,j] <= 4;
### C1 ###
s.t. L1argeCircle_for_C1 {i in C1, j in IP}:
((xr[i,j]-100)^2 + (yr[i,j]-100)^2) / 8014.224 >= 1;
s.t. M1iddleCircle_for_C1 {i in C1, j in IP}:
((xr[i,j]-56)^2 + (yr[i,j]-56)^2) /  2072.271 >= 1;
s.t. S1mallCircle_for_C1 {i in C1, j in IP}:
((xr[i,j]-34)^2 + (yr[i,j]-34)^2) / 553.294 >= 1;
s.t. L2argeCircle_for_C1 {i in C1, j in IP}:
((xf[i,j]-100)^2 + (yf[i,j]-100)^2) / 8014.224 >= 1;
s.t. M2iddleCircle_for_C1 {i in C1, j in IP}:
((xf[i,j]-56)^2 + (yf[i,j]-56)^2) /  2072.271 >= 1;
s.t. S2mallCircle_for_C1 {i in C1, j in IP}:
((xf[i,j]-34)^2 + (yf[i,j]-34)^2) / 553.294 >= 1;
s.t. limit_xr_C1 {i in C1, j in IP}:
xr[i,j] - R >= -12;
s.t. limit_xf_C1 {i in C1, j in IP}:
xf[i,j] - R >= -12;
s.t. limit_yr_C1 {i in C1, j in IP}:
yr[i,j] - R >= -12;
s.t. limit_yf_C1 {i in C1, j in IP}:
yf[i,j] - R >= -12;
s.t. terminal_x_C1 {m in C1}:
xr[m,Nfe] - R >= 12;
s.t. terminal_xr_C1 {m in C1}:
yr[m,Nfe] + R <= 0;
s.t. terminal_xf_C1 {m in C1}:
yf[m,Nfe] + R <= 0;
s.t. interval_theta_C1 {i in C1, j in IPmnt}:
-2 <= theta[i,j] <= 1;
### C2 ###
s.t. limit_xr_C2 {i in C2, j in IP}:
-12 + R <= xr[i,j] <= -R;
s.t. limit_xf_C2 {i in C2, j in IP}:
-12 + R <= xf[i,j] <= -R;
s.t. terminal_x_C2 {m in C2}:
yr[m,Nfe] + R <= -12;
s.t. interval_theta_C2 {i in C2, j in IPmnt}:
-2 <= theta[i,j] <= -1;
### C3 ###
s.t. L1argeCircle_for_C3 {i in C3, j in IP}:
((xr[i,j]+100)^2 + (yr[i,j]-100)^2) / 8014.224 >= 1;
s.t. M1iddleCircle_for_C3 {i in C3, j in IP}:
((xr[i,j]+56)^2 + (yr[i,j]-56)^2) /  2072.271 >= 1;
s.t. S1mallCircle_for_C3 {i in C3, j in IP}:
((xr[i,j]+34)^2 + (yr[i,j]-34)^2) / 553.294 >= 1;
s.t. L2argeCircle_for_C3 {i in C3, j in IP}:
((xf[i,j]+100)^2 + (yf[i,j]-100)^2) / 8014.224 >= 1;
s.t. M2iddleCircle_for_C3 {i in C3, j in IP}:
((xf[i,j]+56)^2 + (yf[i,j]-56)^2) /  2072.271 >= 1;
s.t. S2mallCircle_for_C3 {i in C3, j in IP}:
((xf[i,j]+34)^2 + (yf[i,j]-34)^2) / 553.294 >= 1;
s.t. limit_xr_C3 {i in C3, j in IP}:
xr[i,j] + R <= 0;
s.t. limit_xf_C3 {i in C3, j in IP}:
xf[i,j] + R <= 0;
s.t. limit_yr_C3 {i in C3, j in IP}:
yr[i,j] - R >= 0;
s.t. limit_yf_C3 {i in C3, j in IP}:
yf[i,j] - R >= 0;
s.t. terminal_x_C3 {m in C3}:
xr[m,Nfe] + R <= -12;
s.t. interval_theta_C3 {i in C3, j in IPmnt}:
-4 <= theta[i,j] <= -1;

s.t. VehicleItoJfr {i in {1..(Nv-1)}, j in {(i+1)..Nv}, kk in I}:
(xf[i,kk] - xr[j,kk])^2 + (yf[i,kk] - yr[j,kk])^2 >= FourSqR;
s.t. VehicleItoJrf {i in {1..(Nv-1)}, j in {(i+1)..Nv}, kk in I}:
(xr[i,kk] - xf[j,kk])^2 + (yr[i,kk] - yf[j,kk])^2 >= FourSqR;
s.t. VehicleItoJff {i in {1..(Nv-1)}, j in {(i+1)..Nv}, kk in I}:
(xf[i,kk] - xf[j,kk])^2 + (yf[i,kk] - yf[j,kk])^2 >= FourSqR;
s.t. VehicleItoJrr {i in {1..(Nv-1)}, j in {(i+1)..Nv}, kk in I}:
(xr[i,kk] - xr[j,kk])^2 + (yr[i,kk] - yr[j,kk])^2 >= FourSqR;

s.t. EQ_init_x {m in {1..Nv}}:
x[m,1] = BV[m,1];
s.t. EQ_init_y {m in {1..Nv}}:
y[m,1] = BV[m,2];
s.t. EQ_init_theta {m in {1..Nv}}:
theta[m,1] = BV[m,3];
s.t. EQ_init_v {m in {1..Nv}}:
v[m,1] = vcommon;
s.t. EQ_init_a {m in {1..Nv}}:
a[m,1] = 0;
s.t. EQ_init_phy {m in {1..Nv}}:
phy[m,1] = 0;
s.t. EQ_init_w {m in {1..Nv}}:
w[m,1] = 0;
s.t. EQ_end_theta {m in {1..Nv}}:
theta[m,Nfe] = BV[m,6];
s.t. EQ_end_v {m in {1..Nv}}:
v[m,Nfe] = vcommon;
s.t. EQ_end_a {m in {1..Nv}}:
a[m,Nfe] = 0;
s.t. EQ_end_phy {m in {1..Nv}}:
phy[m,Nfe] = 0;
s.t. EQ_end_w {m in {1..Nv}}:
w[m,Nfe] = 0;

data;
param: cav_index_start := include cav_index_start;
param: cav_index_end := include cav_index_end;
param: BV := include BV;
param: BP := include BP;
param: TimeLine := include TimeLine;