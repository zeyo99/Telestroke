# Parameters
param m;          # Number of Hub Hospitals
param n;          # Number of Spoke Hospitals

#set
set M := 1..m;
set N := 1..n;

#Government參數
param alpha;
param beta;
param S;

#hub參數
#param s {1..m};    變成gov變數
param gamma;            # Cost per patient for each HUB hospital
param theta1;			   # Predefined ratio for payment-cost constraint
param theta2;
param c_hub {1..m};        # Additional cost for each HUB Hospital
param B_hub {1..m};        # Budget for each HUB Hospital

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

#Gov變數
var s {i in M} >= 0;

#Hub變數 
var p {i in M} >= 0;
var d {i in M} >= 0;

#Spoke變數 (所有 hub 共用，只寫一次)
var x {i in 1..m, j in 1..n} >= 0;
var Y_p{j in 1..n} >=0;
var Y_n{j in 1..n} >=0;
var Z_p{i in 1..m, j in 1..n-1, k in j+1..n} >= 0;
var Z_n{i in 1..m, j in 1..n-1, k in j+1..n} >= 0;

#Spoke對偶變數(拉格朗日乘數)
# 上下界由 Spoke stationarity 推得 (有效界限，不會切掉可行解)
var lam{i in 1..m} >= 0, <= M_penalty + (n-1)*delta;
var mu{j in  1..n} >= 0, <= (M_penalty + (n-1)*delta) / c_spoke[j];
var L{j in 1..n} >= -M_penalty, <= M_penalty;
var R{i in 1..m, j in 1..n-1, k in j+1..n} >= -delta, <= delta;

var xd{i in 1..m, j in 1..n} >= 0;
var Yd_p{j in 1..n} >= 0, <= 2*M_penalty;
var Yd_n{j in 1..n} >= 0, <= 2*M_penalty;
var Zd_p{i in 1..m, j in 1..n-1, k in j+1..n} >= 0, <= 2*delta;
var Zd_n{i in 1..m, j in 1..n-1, k in j+1..n} >= 0, <= 2*delta;

#Spoke slack變數
var S_c{i in M} >= 0;
var S_b{j in N} >= 0;

#--- 各 Hub 之 Lagrange 乘數 (第一個索引 h = 哪個 Hub 的 MPEC) ---
# Hub 限制式 (>= 0)
var rho   {h in M} >= 0;        # budget_constraint_hub[h]
var sigma {h in M} >= 0;        # ratio_constraints_hub[h]
var eta_d {h in M} >= 0;        # d[h] >= 0
var eta_p {h in M} >= 0;        # p[h] >= 0

# Spoke Stationarity 之乘數 W (free)
var W_x  {h in M, i in 1..m, j in 1..n};
var W_yp {h in M, j in 1..n};
var W_yn {h in M, j in 1..n};
var W_zp {h in M, i in 1..m, j in 1..n-1, k in j+1..n};
var W_zn {h in M, i in 1..m, j in 1..n-1, k in j+1..n};

# Spoke 原限制式之乘數 zeta (free)
var zeta_c {h in M, i in 1..m};
var zeta_B {h in M, j in 1..n};
var zeta_y {h in M, j in 1..n};
var zeta_z {h in M, i in 1..m, j in 1..n-1, k in j+1..n};

# Spoke 變數非負限制之乘數 V (free, M-stationarity)
var V_x   {h in M, i in 1..m, j in 1..n};
var V_dx  {h in M, i in 1..m, j in 1..n};
var V_yp  {h in M, j in 1..n};
var V_dyp {h in M, j in 1..n};
var V_yn  {h in M, j in 1..n};
var V_dyn {h in M, j in 1..n};
var V_zp  {h in M, i in 1..m, j in 1..n-1, k in j+1..n};
var V_dzp {h in M, i in 1..m, j in 1..n-1, k in j+1..n};
var V_zn  {h in M, i in 1..m, j in 1..n-1, k in j+1..n};
var V_dzn {h in M, i in 1..m, j in 1..n-1, k in j+1..n};
var V_lam {h in M, i in 1..m};
var V_sc  {h in M, i in 1..m};
var V_mu  {h in M, j in 1..n};
var V_sB  {h in M, j in 1..n};

#--- 互補用二元變數 (b = 1 代表該變數 = 0) ---
# Hub level (各 Hub 一組)
var b_rho   {h in M} binary;     # 1 => rho[h] = 0
var b_sigma {h in M} binary;     # 1 => sigma[h] = 0
var b_etad  {h in M} binary;     # 1 => eta_d[h] = 0
var b_etap  {h in M} binary;     # 1 => eta_p[h] = 0

# Spoke level (Spoke 的 active set 只有一個，所有 Hub 共用)
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
# 目標式：醫療覆蓋率最大化
#==================================================================
maximize GOV:
    alpha * sum{i in M} sum{j in N} x[i,j] 
    - beta * sum{i in M} s[i] * sum{j in N} x[i,j] ;

#==================================================================
# 限制式：政府預算
#==================================================================

s.t. gov_budget:
    sum{i in M} s[i] * sum{j in N} x[i,j] <= S;

#=======================TwoHub_MPEC_kkt============================

#==================================================================
# ① Stationarity (各 Hub 之 MPEC 各一組)
#==================================================================

#--- (Hub 變數) ---
s.t. Stat_d {h in M}:
    - c_hub[h] - rho[h] * c_hub[h] + theta2 * sigma[h] + eta_d[h] + zeta_c[h,h] = 0;

s.t. Stat_p {h in M}:
    sum {j in 1..n} x[h,j] + sigma[h] + eta_p[h]
    - sum {j in 1..n} W_x[h,h,j] * (mu[j] + 1)
    - sum {j in 1..n} zeta_B[h,j] * x[h,j] = 0;

#--- (Spoke 原變數) ---
# i = h 有收入項 (p[h] + s[h] + gamma*s[h]*rho[h])；i != h 無收入項，僅透過 zeta 影響 Spoke 可行域
s.t. Stat_x {h in M, i in 1..m, j in 1..n}:
    (if i == h then (p[h] + s[h] + gamma * s[h] * rho[h]) else 0)
    - zeta_c[h,i]
    - (p[i] + c_spoke[j]) * zeta_B[h,j]
    - zeta_y[h,j]
    - sum {k in j+1..n} zeta_z[h,i,j,k]
    + sum {k in 1..j-1} zeta_z[h,i,k,j]
    + V_x[h,i,j] = 0;

#這裡都有index h 每個hub都對原FOC作互補差額 手寫版要補上
s.t. Stat_Yp {h in M, j in 1..n}: - zeta_y[h,j] + V_yp[h,j] = 0;
s.t. Stat_Yn {h in M, j in 1..n}:   zeta_y[h,j] + V_yn[h,j] = 0;

s.t. Stat_Zp {h in M, i in 1..m, j in 1..n-1, k in j+1..n}: - zeta_z[h,i,j,k] + V_zp[h,i,j,k] = 0;
s.t. Stat_Zn {h in M, i in 1..m, j in 1..n-1, k in j+1..n}:   zeta_z[h,i,j,k] + V_zn[h,i,j,k] = 0;

s.t. Stat_Sc {h in M, i in 1..m}: - zeta_c[h,i] + V_sc[h,i] = 0;
s.t. Stat_Sb {h in M, j in 1..n}: - zeta_B[h,j] + V_sB[h,j] = 0;
                                                
#--- (Spoke 對偶變數) ---
s.t. Stat_lam {h in M, i in 1..m}:
    - sum {j in 1..n} W_x[h,i,j] + V_lam[h,i] = 0;

s.t. Stat_mu {h in M, j in 1..n}:
    - sum {i in 1..m} W_x[h,i,j] * (p[i] + c_spoke[j]) + V_mu[h,j] = 0;
                                    #是i還是h 待確認
s.t. Stat_L {h in M, j in 1..n}:
    - sum {i in 1..m} W_x[h,i,j] - W_yp[h,j] + W_yn[h,j] = 0;

s.t. Stat_R {h in M, i in 1..m, j in 1..n-1, k in j+1..n}:
    - W_x[h,i,j] + W_x[h,i,k] - W_zp[h,i,j,k] + W_zn[h,i,j,k] = 0;

s.t. Stat_xd  {h in M, i in 1..m, j in 1..n}: W_x[h,i,j] + V_dx[h,i,j] = 0;
s.t. Stat_Ydp {h in M, j in 1..n}: W_yp[h,j] + V_dyp[h,j] = 0;
s.t. Stat_Ydn {h in M, j in 1..n}: W_yn[h,j] + V_dyn[h,j] = 0;
s.t. Stat_Zdp {h in M, i in 1..m, j in 1..n-1, k in j+1..n}: W_zp[h,i,j,k] + V_dzp[h,i,j,k] = 0;
s.t. Stat_Zdn {h in M, i in 1..m, j in 1..n-1, k in j+1..n}: W_zn[h,i,j,k] + V_dzn[h,i,j,k] = 0;

#==================================================================
# ② Complementary Slackness
#==================================================================

#--- (Hub 限制式，各 Hub 一組) ---
# 原限制式已保證 slack >= 0，故 b=0 時只需加反向不等式即可迫使 slack = 0
# (避免 "==> 等式" 被 AMPL 拆成 or 而產生 aux var tolerance violation)
# rho[h] ⊥ (B_hub + gamma*s*sum x[h,j] - c_hub*d[h])
s.t. Comp_rho_1 {h in M}: b_rho[h] == 1 ==> rho[h] = 0;
s.t. Comp_rho_2 {h in M}: b_rho[h] == 0 ==> c_hub[h] * d[h] >= B_hub[h] + gamma * s[h] * sum {j in 1..n} x[h,j];

# sigma[h] ⊥ (p[h] + theta2*d[h] - theta1*c_hub)
s.t. Comp_sigma_1 {h in M}: b_sigma[h] == 1 ==> sigma[h] = 0;
s.t. Comp_sigma_2 {h in M}: b_sigma[h] == 0 ==> p[h] + theta2 * d[h] <= theta1 * c_hub[h];

# eta_p[h] ⊥ p[h] ; eta_d[h] ⊥ d[h]
s.t. Comp_etap_1 {h in M}: b_etap[h] == 1 ==> eta_p[h] = 0;
s.t. Comp_etap_2 {h in M}: b_etap[h] == 0 ==> p[h] <= 0;
s.t. Comp_etad_1 {h in M}: b_etad[h] == 1 ==> eta_d[h] = 0;
s.t. Comp_etad_2 {h in M}: b_etad[h] == 0 ==> d[h] <= 0;

#--- (Spoke level 變數互補：寫一次 + 各 Hub 的 V 互補, M-stationarity) ---
# x ⊥ xd
#x ⊥ Dx
s.t. Cx_1 {i in 1..m, j in 1..n}: b_x[i,j]  == 1 ==> x[i,j]  = 0;
s.t. Cx_2 {i in 1..m, j in 1..n}: b_xd[i,j] == 1 ==> xd[i,j] = 0;
s.t. Cx_3 {i in 1..m, j in 1..n}: b_x[i,j] + b_xd[i,j] >= 1;
# Vx·x = 0、Vdx·Dx = 0
s.t. Cx_4 {h in M, i in 1..m, j in 1..n}: b_x[i,j]  == 0 ==> V_x[h,i,j]  = 0;
s.t. Cx_5 {h in M, i in 1..m, j in 1..n}: b_xd[i,j] == 0 ==> V_dx[h,i,j] = 0;
#biactive : x= 0, Dx=0 => Vx >= 0, Vdx >= 0
s.t. Cx_6 {h in M, i in 1..m, j in 1..n}: b_xd[i,j] == 1 ==> V_x[h,i,j]  >= 0;
s.t. Cx_7 {h in M, i in 1..m, j in 1..n}: b_x[i,j]  == 1 ==> V_dx[h,i,j] >= 0;

# Y_p ⊥ Yd_p
s.t. Cyp_1 {j in 1..n}: b_yp[j]  == 1 ==> Y_p[j]  = 0;
s.t. Cyp_2 {j in 1..n}: b_ydp[j] == 1 ==> Yd_p[j] = 0;
s.t. Cyp_3 {j in 1..n}: b_yp[j] + b_ydp[j] >= 1;
s.t. Cyp_4 {h in M, j in 1..n}: b_yp[j]  == 0 ==> V_yp[h,j]  = 0;
s.t. Cyp_5 {h in M, j in 1..n}: b_ydp[j] == 0 ==> V_dyp[h,j] = 0;
s.t. Cyp_6 {h in M, j in 1..n}: b_ydp[j] == 1 ==> V_yp[h,j]  >= 0;
s.t. Cyp_7 {h in M, j in 1..n}: b_yp[j]  == 1 ==> V_dyp[h,j] >= 0;

# Y_n ⊥ Yd_n
s.t. Cyn_1 {j in 1..n}: b_yn[j]  == 1 ==> Y_n[j]  = 0;
s.t. Cyn_2 {j in 1..n}: b_ydn[j] == 1 ==> Yd_n[j] = 0;
s.t. Cyn_3 {j in 1..n}: b_yn[j] + b_ydn[j] >= 1;
s.t. Cyn_4 {h in M, j in 1..n}: b_yn[j]  == 0 ==> V_yn[h,j]  = 0;
s.t. Cyn_5 {h in M, j in 1..n}: b_ydn[j] == 0 ==> V_dyn[h,j] = 0;
s.t. Cyn_6 {h in M, j in 1..n}: b_ydn[j] == 1 ==> V_yn[h,j]  >= 0;
s.t. Cyn_7 {h in M, j in 1..n}: b_yn[j]  == 1 ==> V_dyn[h,j] >= 0;

# Z_p ⊥ Zd_p
s.t. Czp_1 {i in 1..m, j in 1..n-1, k in j+1..n}: b_zp[i,j,k]  == 1 ==> Z_p[i,j,k]  = 0;
s.t. Czp_2 {i in 1..m, j in 1..n-1, k in j+1..n}: b_zdp[i,j,k] == 1 ==> Zd_p[i,j,k] = 0;
s.t. Czp_3 {i in 1..m, j in 1..n-1, k in j+1..n}: b_zp[i,j,k] + b_zdp[i,j,k] >= 1;
s.t. Czp_4 {h in M, i in 1..m, j in 1..n-1, k in j+1..n}: b_zp[i,j,k]  == 0 ==> V_zp[h,i,j,k]  = 0;
s.t. Czp_5 {h in M, i in 1..m, j in 1..n-1, k in j+1..n}: b_zdp[i,j,k] == 0 ==> V_dzp[h,i,j,k] = 0;
s.t. Czp_6 {h in M, i in 1..m, j in 1..n-1, k in j+1..n}: b_zdp[i,j,k] == 1 ==> V_zp[h,i,j,k]  >= 0;
s.t. Czp_7 {h in M, i in 1..m, j in 1..n-1, k in j+1..n}: b_zp[i,j,k]  == 1 ==> V_dzp[h,i,j,k] >= 0;

# Z_n ⊥ Zd_n
s.t. Czn_1 {i in 1..m, j in 1..n-1, k in j+1..n}: b_zn[i,j,k]  == 1 ==> Z_n[i,j,k]  = 0;
s.t. Czn_2 {i in 1..m, j in 1..n-1, k in j+1..n}: b_zdn[i,j,k] == 1 ==> Zd_n[i,j,k] = 0;
s.t. Czn_3 {i in 1..m, j in 1..n-1, k in j+1..n}: b_zn[i,j,k] + b_zdn[i,j,k] >= 1;
s.t. Czn_4 {h in M, i in 1..m, j in 1..n-1, k in j+1..n}: b_zn[i,j,k]  == 0 ==> V_zn[h,i,j,k]  = 0;
s.t. Czn_5 {h in M, i in 1..m, j in 1..n-1, k in j+1..n}: b_zdn[i,j,k] == 0 ==> V_dzn[h,i,j,k] = 0;
s.t. Czn_6 {h in M, i in 1..m, j in 1..n-1, k in j+1..n}: b_zdn[i,j,k] == 1 ==> V_zn[h,i,j,k]  >= 0;
s.t. Czn_7 {h in M, i in 1..m, j in 1..n-1, k in j+1..n}: b_zn[i,j,k]  == 1 ==> V_dzn[h,i,j,k] >= 0;

# lam ⊥ S_c
s.t. Clam_1 {i in 1..m}: b_lam[i] == 1 ==> lam[i] = 0;
s.t. Clam_2 {i in 1..m}: b_sc[i]  == 1 ==> S_c[i] = 0;
s.t. Clam_3 {i in 1..m}: b_lam[i] + b_sc[i] >= 1;
s.t. Clam_4 {h in M, i in 1..m}: b_lam[i] == 0 ==> V_lam[h,i] = 0;
s.t. Clam_5 {h in M, i in 1..m}: b_sc[i]  == 0 ==> V_sc[h,i]  = 0;
s.t. Clam_6 {h in M, i in 1..m}: b_sc[i]  == 1 ==> V_lam[h,i] >= 0;
s.t. Clam_7 {h in M, i in 1..m}: b_lam[i] == 1 ==> V_sc[h,i]  >= 0;

# mu ⊥ S_b
s.t. Cmu_1 {j in 1..n}: b_mu[j] == 1 ==> mu[j]  = 0;
s.t. Cmu_2 {j in 1..n}: b_sb[j] == 1 ==> S_b[j] = 0;
s.t. Cmu_3 {j in 1..n}: b_mu[j] + b_sb[j] >= 1;
s.t. Cmu_4 {h in M, j in 1..n}: b_mu[j] == 0 ==> V_mu[h,j] = 0;
s.t. Cmu_5 {h in M, j in 1..n}: b_sb[j] == 0 ==> V_sB[h,j] = 0;
s.t. Cmu_6 {h in M, j in 1..n}: b_sb[j] == 1 ==> V_mu[h,j] >= 0;
s.t. Cmu_7 {h in M, j in 1..n}: b_mu[j] == 1 ==> V_sB[h,j] >= 0;

#==================================================================
# ③ 原限制式
#==================================================================

#--- A. Hub 限制式 (各 Hub 一組) ---
subject to budget_constraint_hub {h in M}:
     c_hub[h] * d[h] <= B_hub[h] + gamma * s[h] * sum{j in 1..n} x[h,j];

subject to ratio_constraints_hub {h in M}:
	p[h] + theta2 * d[h] >= theta1 * c_hub[h];

#--- B-i. Stationarity of Spoke (共用，寫一次) ---
s.t. Stationarity_x {i in 1..m, j in 1..n}:
     lam[i] + (mu[j]+ 1) * (p[i] + c_spoke[j]) + pref_noise[i,j] + L[j]
     + sum {k in j+1..n} R[i,j,k] - sum {k in 1..j-1} R[i,k,j] - xd[i,j] = 0;

s.t. StationarityForY_p{j in 1..n}:
	M_penalty + L[j] - Yd_p[j] = 0;

s.t. StationarityForY_n{j in 1..n}:
	M_penalty - L[j] - Yd_n[j] = 0;

s.t. StationarityForZ_p{i in 1..m, j in 1..n-1, k in j+1..n}:
	delta + R[i,j,k] - Zd_p[i,j,k] = 0;

s.t. StationarityForZ_n{i in 1..m, j in 1..n-1, k in j+1..n}:
	delta - R[i,j,k] - Zd_n[i,j,k] = 0;

#--- B-ii. Spoke 之原限制式 (共用，寫一次) ---
s.t. capacity_constraint_slack {i in 1..m}:
    sum {j in 1..n} x[i,j] + S_c[i] = d[i];

s.t. cost_constraint_slack {j in 1..n}:
    sum {i in 1..m} (p[i] + c_spoke[j]) * x[i,j] + S_b[j] = B_spoke[j];

s.t. y_var_constraint{j in 1..n}:
	Y_p[j] - Y_n[j] = T[j] - sum{i in 1..m}x[i,j];

s.t. z_var_constraint {i in 1..m, j in 1..n-1, k in j+1..n}:
    Z_p[i,j,k] - Z_n[i,j,k] = x[i,k] - x[i,j] ;

# ④ 變數本身之限制：已寫在 var 宣告 (>= 0 / free)
