# Number of Hub Hospitals
param m;

# Number of spoke Hospitals
param n;

#set
set M := 1..m;
set N := 1..n;

# Demand of spoke Hospitals T_j
param T {1..n};

# Capacity of HUB Hospitals d_i
param d {1..m};

# Cost per unit of capacity for each HUB Hospital p_i
param p {1..m};

# Additional cost for each spoke Hospital c_j
param c_spoke {1..n};

# Budget limit for each spoke Hospital B_j
param B_spoke {1..n};

# Regularization parameter delta
param delta;

#大M懲罰值(防止trivial solution)
param M_penalty;

#異質性偏好(多個hub價格相同時spoke能實現隨機選擇)
param pref_noise {1..m, 1..n} ; #default Uniform(0.1, 5.0);

#Variables
var x {i in 1..m, j in 1..n} >= 0; #convex programming 要求變數連續可微
var Y_p{j in 1..n} >=0;
var Y_n{j in 1..n} >=0; 
var Z_p{i in 1..m, j in 1..n-1, k in j+1..n} >= 0;
var Z_n{i in 1..m, j in 1..n-1, k in j+1..n} >= 0;

# dual vars
var lam{i in 1..m} >=0;
var mu{j in  1..n} >=0;
var L{j in 1..n} ;
var R{i in 1..m, j in 1..n-1, k in j+1..n} ;

var xd{i in 1..m, j in 1..n} >= 0;
var Yd_p{j in 1..n} >=0;
var Yd_n{j in 1..n} >=0; 
var Zd_p{i in 1..m, j in 1..n-1, k in j+1..n} >= 0;
var Zd_n{i in 1..m, j in 1..n-1, k in j+1..n} >= 0;

#slack變數
var S_c{i in M} >= 0;
var S_b{j in N} >= 0;

#互補二元變數
var DV_cap    {1..m} binary;
var DV_budget {1..n} binary;
var DV_x      {1..m, 1..n} binary;
var DV_yp     {1..n} binary;
var DV_yn     {1..n} binary;
var DV_zp     {i in 1..m,j in 1..n-1, k in j+1..n} binary; # 注意 k 的範圍跟隨 j
var DV_zn     {i in 1..m,j in 1..n-1, k in j+1..n} binary;

# Objective Function 
minimize total_obj:
       M_penalty * sum {j in 1..n}(Y_p[j] + Y_n[j]) 
       + delta * ( sum{i in 1..m} sum{j in 1..n-1} sum{k in j+1..n} (Z_p[i,j,k] + Z_n[i,j,k]))
       + sum{j in N} sum{i in M} ( (p[i] + c_spoke[j] + pref_noise[i,j] ) * x[i,j] ); 
                              

# KKT Stationarity(FOC 條件) 
# ==========================================
s.t. Stationarity_x {i in 1..m, j in 1..n}:
     lam[i]+ (mu[j]+1) * (p[i] + c_spoke[j]) + pref_noise[i,j] + L[j] 
      + sum {k in j+1..n} R[i,j,k] - sum {k in 1..j-1} R[i,k,j] - xd[i,j] = 0; 

s.t. StationarityForY_p{j in 1..n}:
	M_penalty + L[j] - Yd_p[j] = 0;
    
s.t. StationarityForY_n{j in 1..n}:
	M_penalty - L[j] - Yd_n[j] = 0;
	
s.t. StationarityForZ_p{i in 1..m, j in 1..n-1, k in j+1..n}:
	delta + R[i,j,k] - Zd_p[i,j,k] = 0;
	
s.t. StationarityForZ_n{i in 1..m, j in 1..n-1, k in j+1..n}:
	delta - R[i,j,k] - Zd_n[i,j,k] = 0;
	
#Comlementary Slackness 替換


 
 #邏輯指示限制式 (Indicator Constraints)
# --- Capacity 互補 (lam * S_c = 0) ---
s.t. Ind_capacity_1 {i in 1..m}: DV_cap[i] == 0 ==> lam[i] == 0;
s.t. Ind_capacity_2 {i in 1..m}: DV_cap[i] == 1 ==> S_c[i] == 0;

# --- Budget 互補 (mu * S_b = 0) ---
s.t. Ind_budget_1 {j in 1..n}: DV_budget[j] == 0 ==> mu[j] == 0;
s.t. Ind_budget_2 {j in 1..n}: DV_budget[j] == 1 ==> S_b[j] == 0;

# --- x 變數互補 (x * xd = 0) ---
s.t. Ind_x_1 {i in 1..m, j in 1..n}: DV_x[i,j] == 0 ==> x[i,j] == 0;
s.t. Ind_x_2 {i in 1..m, j in 1..n}: DV_x[i,j] == 1 ==> xd[i,j] == 0;

# --- Y_p 變數互補 (Y_p * Yd_p = 0) ---
s.t. Ind_yp_1 {j in 1..n}: DV_yp[j] == 0 ==> Y_p[j] == 0;
s.t. Ind_yp_2 {j in 1..n}: DV_yp[j] == 1 ==> Yd_p[j] == 0;

# --- Y_n 變數互補 (Y_n * Yd_n = 0) ---
s.t. Ind_yn_1 {j in 1..n}: DV_yn[j] == 0 ==> Y_n[j] == 0;
s.t. Ind_yn_2 {j in 1..n}: DV_yn[j] == 1 ==> Yd_n[j] == 0;

# --- Z_p 變數互補 (Z_p * Zd_p = 0) ---
s.t. Ind_zp_1 {i in 1..m, j in 1..n-1, k in j+1..n}: DV_zp[i,j,k] == 0 ==> Z_p[i,j,k] == 0;
s.t. Ind_zp_2 {i in 1..m, j in 1..n-1, k in j+1..n}: DV_zp[i,j,k] == 1 ==> Zd_p[i,j,k] == 0;

# --- Z_n 變數互補 (Z_n * Zd_n = 0) ---
s.t. Ind_zn_1 {i in 1..m, j in 1..n-1, k in j+1..n}: DV_zn[i,j,k] == 0 ==> Z_n[i,j,k] == 0;
s.t. Ind_zn_2 {i in 1..m, j in 1..n-1, k in j+1..n}: DV_zn[i,j,k] == 1 ==> Zd_n[i,j,k] == 0;

# 原限制式 
# ==========================================


s.t. capacity_constraint {i in 1..m}:
    sum {j in 1..n} x[i,j] + S_c[i] = d[i];   
 

s.t. cost_constraint_slack {j in 1..n}:
    sum {i in 1..m} (p[i] + c_spoke[j]) * x[i,j] + S_b[j] = B_spoke[j];
    
s.t. y_var_constraint{j in 1..n}:
	Y_p[j] - Y_n[j] = T[j] - sum{i in 1..m}x[i,j];


s.t. z_var_constraint {i in 1..m, j in 1..n-1, k in j+1..n}:
    Z_p[i,j,k] - Z_n[i,j,k] = x[i,k] - x[i,j] ; 
    