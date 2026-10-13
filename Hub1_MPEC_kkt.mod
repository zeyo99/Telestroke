# Hub 1 MPEC 的 KKT 系統 (對應 Hub_mpec_kkt_2.1.pdf)
# Hub 2 的 p2, d2 視為參數

# Parameters
param m;          # Number of Hub Hospitals
param n;          # Number of Spoke Hospitals

#set
set M := 1..m;
set N := 1..n;

#hub參數
param s {1..m};
param gamma;            # Cost per patient for each HUB hospital
param theta1;			   # Predefined ratio for payment-cost constraint
param theta2;
param c_hub {1..m};        # Additional cost for each HUB Hospital
param B_hub {1..m};        # Budget for each HUB Hospital

#固定hub2的變數為參數
param p2;
param d2;

#spoke參數
param c_spoke {1..n};      # Additional cost for each spoke Hospital
param B_spoke {1..n};      # Budget for each spoke Hospital
param T {1..n};            # Demand of Spoke Hospitals
param delta;               # Regularization parameter

#大M懲罰值(防止trivial solution)
param M_penalty;

#異質性偏好(多個hub價格相同時spoke能實現隨機選擇)
param pref_noise {1..m, 1..n} ; #default Uniform(0.1, 5.0);

#==================================================================
# Variables
#==================================================================

#Hub 1變數
var p1 >= 0;
var d1 >= 0;

#Spoke變數
var x {i in 1..m, j in 1..n} >= 0;
var Y_p{j in 1..n} >=0;
var Y_n{j in 1..n} >=0;
var Z_p{i in 1..m, j in 1..n-1, k in j+1..n} >= 0;
var Z_n{i in 1..m, j in 1..n-1, k in j+1..n} >= 0;

#Spoke對偶變數(拉格朗日乘數)
var lam{i in 1..m} >=0;
var mu{j in  1..n} >=0;
var L{j in 1..n};
var R{i in 1..m, j in 1..n-1, k in j+1..n};

var xd{i in 1..m, j in 1..n} >= 0;
var Yd_p{j in 1..n} >=0;
var Yd_n{j in 1..n} >=0;
var Zd_p{i in 1..m, j in 1..n-1, k in j+1..n} >= 0;
var Zd_n{i in 1..m, j in 1..n-1, k in j+1..n} >= 0;

#Spoke slack變數
var S_c{i in M} >= 0;
var S_b{j in N} >= 0;

#--- Hub 1 Lagrange 乘數 ---
# Hub 限制式 (>= 0)
var rho1   >= 0;        # budget_constraint_hub
var sigma1 >= 0;        # ratio_constraints_hub
var eta_d  >= 0;        # d1 >= 0
var eta_p  >= 0;        # p1 >= 0

# Spoke Stationarity 之乘數 W (free)
var W_x  {i in 1..m, j in 1..n};
var W_yp {j in 1..n};
var W_yn {j in 1..n};
var W_zp {i in 1..m, j in 1..n-1, k in j+1..n};
var W_zn {i in 1..m, j in 1..n-1, k in j+1..n};

# Spoke 原限制式之乘數 zeta (free)
var zeta_c {i in 1..m};
var zeta_B {j in 1..n};
var zeta_y {j in 1..n};
var zeta_z {i in 1..m, j in 1..n-1, k in j+1..n};

# Spoke 變數非負限制之乘數 V (free, M-stationarity)
var V_x   {i in 1..m, j in 1..n};
var V_dx  {i in 1..m, j in 1..n};
var V_yp  {j in 1..n};
var V_dyp {j in 1..n};
var V_yn  {j in 1..n};
var V_dyn {j in 1..n};
var V_zp  {i in 1..m, j in 1..n-1, k in j+1..n};
var V_dzp {i in 1..m, j in 1..n-1, k in j+1..n};
var V_zn  {i in 1..m, j in 1..n-1, k in j+1..n};
var V_dzn {i in 1..m, j in 1..n-1, k in j+1..n};
var V_lam {i in 1..m};
var V_sc  {i in 1..m};
var V_mu  {j in 1..n};
var V_sB  {j in 1..n};

#--- 互補用二元變數 (b = 1 代表該變數 = 0) ---
# Hub level
var b_rho   binary;     # 1 => rho1 = 0
var b_sigma binary;     # 1 => sigma1 = 0
var b_etad  binary;     # 1 => eta_d = 0
var b_etap  binary;     # 1 => eta_p = 0

# Spoke level
var b_x   {i in 1..m, j in 1..n} binary;
var b_xd  {i in 1..m, j in 1..n} binary;
var b_yp  {j in 1..n} binary;
var b_ydp {j in 1..n} binary;
var b_yn  {j in 1..n} binary;
var b_ydn {j in 1..n} binary;
var b_zp  {i in 1..m, j in 1..n-1, k in j+1..n} binary;
var b_zdp {i in 1..m, j in 1..n-1, k in j+1..n} binary;
var b_zn  {i in 1..m, j in 1..n-1, k in j+1..n} binary;
var b_zdn {i in 1..m, j in 1..n-1, k in j+1..n} binary;
var b_lam {i in 1..m} binary;
var b_sc  {i in 1..m} binary;
var b_mu  {j in 1..n} binary;
var b_sb  {j in 1..n} binary;

#==================================================================
# hub1目標式
#==================================================================
maximize Hub1_profit:
    (p1 + s[1]) * (sum {j in 1..n} x[1,j]) - (c_hub[1] * d1);

#==================================================================
# ① Stationarity
#==================================================================

#--- (Hub 變數) ---
s.t. Stat_d1:
    - c_hub[1] - rho1 * c_hub[1] + theta2 * sigma1 + eta_d + zeta_c[1] = 0;

s.t. Stat_p1:
    sum {j in 1..n} x[1,j] + sigma1 + eta_p
    - sum {j in 1..n} W_x[1,j] * (mu[j] + 1)
    - sum {j in 1..n} zeta_B[j] * x[1,j] = 0;

#--- (Spoke 原變數) ---
# h = 1 有收入項 (p1 + s1 + gamma*s1*rho1)；h = 2 無收入項，僅透過 zeta 影響 Spoke 可行域
s.t. Stat_x {i in 1..m, j in 1..n}:
    (if i == 1 then (p1 + s[1] + gamma * s[1] * rho1) else 0)
    - zeta_c[i]
    - ((if i == 1 then p1 else p2) + c_spoke[j]) * zeta_B[j]
    - zeta_y[j]
    - sum {k in j+1..n} zeta_z[i,j,k]
    + sum {k in 1..j-1} zeta_z[i,k,j]
    + V_x[i,j] = 0;

s.t. Stat_Yp {j in 1..n}: - zeta_y[j] + V_yp[j] = 0;
s.t. Stat_Yn {j in 1..n}:   zeta_y[j] + V_yn[j] = 0;

s.t. Stat_Zp {i in 1..m, j in 1..n-1, k in j+1..n}: - zeta_z[i,j,k] + V_zp[i,j,k] = 0;
s.t. Stat_Zn {i in 1..m, j in 1..n-1, k in j+1..n}:   zeta_z[i,j,k] + V_zn[i,j,k] = 0;

s.t. Stat_Sc {i in 1..m}: - zeta_c[i] + V_sc[i] = 0;
s.t. Stat_Sb {j in 1..n}: - zeta_B[j] + V_sB[j] = 0;

#--- (Spoke 對偶變數) ---
s.t. Stat_lam {i in 1..m}:
    - sum {j in 1..n} W_x[i,j] + V_lam[i] = 0;

s.t. Stat_mu {j in 1..n}:
    - sum {i in 1..m} W_x[i,j] * ((if i == 1 then p1 else p2) + c_spoke[j]) + V_mu[j] = 0;

s.t. Stat_L {j in 1..n}:
    - sum {i in 1..m} W_x[i,j] - W_yp[j] + W_yn[j] = 0;

s.t. Stat_R {i in 1..m, j in 1..n-1, k in j+1..n}:
    - W_x[i,j] + W_x[i,k] - W_zp[i,j,k] + W_zn[i,j,k] = 0;

s.t. Stat_xd  {i in 1..m, j in 1..n}: W_x[i,j] + V_dx[i,j] = 0;
s.t. Stat_Ydp {j in 1..n}: W_yp[j] + V_dyp[j] = 0;
s.t. Stat_Ydn {j in 1..n}: W_yn[j] + V_dyn[j] = 0;
s.t. Stat_Zdp {i in 1..m, j in 1..n-1, k in j+1..n}: W_zp[i,j,k] + V_dzp[i,j,k] = 0;
s.t. Stat_Zdn {i in 1..m, j in 1..n-1, k in j+1..n}: W_zn[i,j,k] + V_dzn[i,j,k] = 0;

#==================================================================
# ② Complementary Slackness
#==================================================================

#--- (Hub 限制式) ---
# 原限制式已保證 slack >= 0，故 b=0 時只需加反向不等式即可迫使 slack = 0
# (避免 "==> 等式" 被 AMPL 拆成 or 而產生 aux var tolerance violation)
# rho1 ⊥ (B_hub + gamma*s1*sum x1j - c_hub*d1)
s.t. Comp_rho_1: b_rho == 1 ==> rho1 = 0;
s.t. Comp_rho_2: b_rho == 0 ==> c_hub[1] * d1 >= B_hub[1] + gamma * s[1] * sum {j in 1..n} x[1,j];

# sigma1 ⊥ (p1 + theta2*d1 - theta1*c_hub)
s.t. Comp_sigma_1: b_sigma == 1 ==> sigma1 = 0;
s.t. Comp_sigma_2: b_sigma == 0 ==> p1 + theta2 * d1 <= theta1 * c_hub[1];

# eta_p ⊥ p1 ; eta_d ⊥ d1
s.t. Comp_etap_1: b_etap == 1 ==> eta_p = 0;
s.t. Comp_etap_2: b_etap == 0 ==> p1 <= 0;
s.t. Comp_etad_1: b_etad == 1 ==> eta_d = 0;
s.t. Comp_etad_2: b_etad == 0 ==> d1 <= 0;

#--- (Spoke level 變數互補 + Hub level V 互補, M-stationarity) ---
# x ⊥ xd
s.t. Cx_1 {i in 1..m, j in 1..n}: b_x[i,j]  == 1 ==> x[i,j]  = 0;
s.t. Cx_2 {i in 1..m, j in 1..n}: b_xd[i,j] == 1 ==> xd[i,j] = 0;
s.t. Cx_3 {i in 1..m, j in 1..n}: b_x[i,j] + b_xd[i,j] >= 1;
s.t. Cx_4 {i in 1..m, j in 1..n}: b_x[i,j]  == 0 ==> V_x[i,j]  = 0;
s.t. Cx_5 {i in 1..m, j in 1..n}: b_xd[i,j] == 0 ==> V_dx[i,j] = 0;
s.t. Cx_6 {i in 1..m, j in 1..n}: b_xd[i,j] == 1 ==> V_x[i,j]  >= 0;
s.t. Cx_7 {i in 1..m, j in 1..n}: b_x[i,j]  == 1 ==> V_dx[i,j] >= 0;

# Y_p ⊥ Yd_p
s.t. Cyp_1 {j in 1..n}: b_yp[j]  == 1 ==> Y_p[j]  = 0;
s.t. Cyp_2 {j in 1..n}: b_ydp[j] == 1 ==> Yd_p[j] = 0;
s.t. Cyp_3 {j in 1..n}: b_yp[j] + b_ydp[j] >= 1;
s.t. Cyp_4 {j in 1..n}: b_yp[j]  == 0 ==> V_yp[j]  = 0;
s.t. Cyp_5 {j in 1..n}: b_ydp[j] == 0 ==> V_dyp[j] = 0;
s.t. Cyp_6 {j in 1..n}: b_ydp[j] == 1 ==> V_yp[j]  >= 0;
s.t. Cyp_7 {j in 1..n}: b_yp[j]  == 1 ==> V_dyp[j] >= 0;

# Y_n ⊥ Yd_n
s.t. Cyn_1 {j in 1..n}: b_yn[j]  == 1 ==> Y_n[j]  = 0;
s.t. Cyn_2 {j in 1..n}: b_ydn[j] == 1 ==> Yd_n[j] = 0;
s.t. Cyn_3 {j in 1..n}: b_yn[j] + b_ydn[j] >= 1;
s.t. Cyn_4 {j in 1..n}: b_yn[j]  == 0 ==> V_yn[j]  = 0;
s.t. Cyn_5 {j in 1..n}: b_ydn[j] == 0 ==> V_dyn[j] = 0;
s.t. Cyn_6 {j in 1..n}: b_ydn[j] == 1 ==> V_yn[j]  >= 0;
s.t. Cyn_7 {j in 1..n}: b_yn[j]  == 1 ==> V_dyn[j] >= 0;

# Z_p ⊥ Zd_p
s.t. Czp_1 {i in 1..m, j in 1..n-1, k in j+1..n}: b_zp[i,j,k]  == 1 ==> Z_p[i,j,k]  = 0;
s.t. Czp_2 {i in 1..m, j in 1..n-1, k in j+1..n}: b_zdp[i,j,k] == 1 ==> Zd_p[i,j,k] = 0;
s.t. Czp_3 {i in 1..m, j in 1..n-1, k in j+1..n}: b_zp[i,j,k] + b_zdp[i,j,k] >= 1;
s.t. Czp_4 {i in 1..m, j in 1..n-1, k in j+1..n}: b_zp[i,j,k]  == 0 ==> V_zp[i,j,k]  = 0;
s.t. Czp_5 {i in 1..m, j in 1..n-1, k in j+1..n}: b_zdp[i,j,k] == 0 ==> V_dzp[i,j,k] = 0;
s.t. Czp_6 {i in 1..m, j in 1..n-1, k in j+1..n}: b_zdp[i,j,k] == 1 ==> V_zp[i,j,k]  >= 0;
s.t. Czp_7 {i in 1..m, j in 1..n-1, k in j+1..n}: b_zp[i,j,k]  == 1 ==> V_dzp[i,j,k] >= 0;

# Z_n ⊥ Zd_n
s.t. Czn_1 {i in 1..m, j in 1..n-1, k in j+1..n}: b_zn[i,j,k]  == 1 ==> Z_n[i,j,k]  = 0;
s.t. Czn_2 {i in 1..m, j in 1..n-1, k in j+1..n}: b_zdn[i,j,k] == 1 ==> Zd_n[i,j,k] = 0;
s.t. Czn_3 {i in 1..m, j in 1..n-1, k in j+1..n}: b_zn[i,j,k] + b_zdn[i,j,k] >= 1;
s.t. Czn_4 {i in 1..m, j in 1..n-1, k in j+1..n}: b_zn[i,j,k]  == 0 ==> V_zn[i,j,k]  = 0;
s.t. Czn_5 {i in 1..m, j in 1..n-1, k in j+1..n}: b_zdn[i,j,k] == 0 ==> V_dzn[i,j,k] = 0;
s.t. Czn_6 {i in 1..m, j in 1..n-1, k in j+1..n}: b_zdn[i,j,k] == 1 ==> V_zn[i,j,k]  >= 0;
s.t. Czn_7 {i in 1..m, j in 1..n-1, k in j+1..n}: b_zn[i,j,k]  == 1 ==> V_dzn[i,j,k] >= 0;

# lam ⊥ S_c
s.t. Clam_1 {i in 1..m}: b_lam[i] == 1 ==> lam[i] = 0;
s.t. Clam_2 {i in 1..m}: b_sc[i]  == 1 ==> S_c[i] = 0;
s.t. Clam_3 {i in 1..m}: b_lam[i] + b_sc[i] >= 1;
s.t. Clam_4 {i in 1..m}: b_lam[i] == 0 ==> V_lam[i] = 0;
s.t. Clam_5 {i in 1..m}: b_sc[i]  == 0 ==> V_sc[i]  = 0;
s.t. Clam_6 {i in 1..m}: b_sc[i]  == 1 ==> V_lam[i] >= 0;
s.t. Clam_7 {i in 1..m}: b_lam[i] == 1 ==> V_sc[i]  >= 0;

# mu ⊥ S_b
s.t. Cmu_1 {j in 1..n}: b_mu[j] == 1 ==> mu[j]  = 0;
s.t. Cmu_2 {j in 1..n}: b_sb[j] == 1 ==> S_b[j] = 0;
s.t. Cmu_3 {j in 1..n}: b_mu[j] + b_sb[j] >= 1;
s.t. Cmu_4 {j in 1..n}: b_mu[j] == 0 ==> V_mu[j] = 0;
s.t. Cmu_5 {j in 1..n}: b_sb[j] == 0 ==> V_sB[j] = 0;
s.t. Cmu_6 {j in 1..n}: b_sb[j] == 1 ==> V_mu[j] >= 0;
s.t. Cmu_7 {j in 1..n}: b_mu[j] == 1 ==> V_sB[j] >= 0;

#==================================================================
# ③ 原限制式
#==================================================================

#--- A. Hub 限制式 ---
subject to budget_constraint_hub :
     c_hub[1] * d1 <= B_hub[1] + gamma * s[1] * sum{j in 1..n} x[1,j];

subject to ratio_constraints_hub :
	p1 + theta2 * d1 >= theta1 * c_hub[1];

#--- B-i. Stationarity of Spoke ---
s.t. Stationarity_x {i in 1..m, j in 1..n}:
     lam[i] + (mu[j]+ 1) * ( (if i==1 then p1 else p2) + c_spoke[j]) + pref_noise[i,j] + L[j]
     + sum {k in j+1..n} R[i,j,k] - sum {k in 1..j-1} R[i,k,j] - xd[i,j] = 0;

s.t. StationarityForY_p{j in 1..n}:
	M_penalty + L[j] - Yd_p[j] = 0;

s.t. StationarityForY_n{j in 1..n}:
	M_penalty - L[j] - Yd_n[j] = 0;

s.t. StationarityForZ_p{i in 1..m, j in 1..n-1, k in j+1..n}:
	delta + R[i,j,k] - Zd_p[i,j,k] = 0;

s.t. StationarityForZ_n{i in 1..m, j in 1..n-1, k in j+1..n}:
	delta - R[i,j,k] - Zd_n[i,j,k] = 0;

#--- B-ii. Spoke 之原限制式 ---
s.t. capacity_constraint_slack {i in 1..m}:
    sum {j in 1..n} x[i,j] + S_c[i] = (if i==1 then d1 else d2);

s.t. cost_constraint_slack {j in 1..n}:
    sum {i in 1..m} ( (if i==1 then p1 else p2) + c_spoke[j]) * x[i,j] + S_b[j] = B_spoke[j];

s.t. y_var_constraint{j in 1..n}:
	Y_p[j] - Y_n[j] = T[j] - sum{i in 1..m}x[i,j];

s.t. z_var_constraint {i in 1..m, j in 1..n-1, k in j+1..n}:
    Z_p[i,j,k] - Z_n[i,j,k] = x[i,k] - x[i,j] ;

# ④ 非負限制：已寫在 var 宣告 (>= 0 / free)
