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

#Variables

#Hub 1變數
var p1 >= 0; 
var d1 >= 0;

#Spoke變數
var x {i in 1..m, j in 1..n} >= 0; #convex programming 要求變數連續可微
var Y_p{j in 1..n} >=0;
var Y_n{j in 1..n} >=0; 
var Z_p{i in 1..m, j in 1..n-1, k in j+1..n} >= 0;
var Z_n{i in 1..m, j in 1..n-1, k in j+1..n} >= 0;

#對偶變數(拉格朗日乘數)
var lam{i in 1..m} >=0;
var mu{j in  1..n} >=0;
var L{j in 1..n};
var R{i in 1..m, j in 1..n-1, k in j+1..n};

var xd{i in 1..m, j in 1..n} >= 0;
var Yd_p{j in 1..n} >=0;
var Yd_n{j in 1..n} >=0; 
var Zd_p{i in 1..m, j in 1..n-1, k in j+1..n} >= 0;
var Zd_n{i in 1..m, j in 1..n-1, k in j+1..n} >= 0;

#slack變數
var S_c{i in M} >= 0;
var S_b{j in N} >= 0;

#邏輯判斷二元變數
var D_cap    {1..m} binary;
var D_budget {1..n} binary;
var D_x      {1..m, 1..n} binary;
var D_yp     {1..n} binary;
var D_zp     {i in 1..m,j in 1..n-1, k in j+1..n} binary; # 注意 k 的範圍跟隨 j
var D_zn     {i in 1..m,j in 1..n-1, k in j+1..n} binary;

#hub1目標式

maximize Hub1_profit:
    (
        # 收入部分：(價格 + 政府補助) * 該醫院收治的所有 Spoke 病患總數
        (p1 +  s[1]) * (sum {j in 1..n} x[1,j]) #雙線性1
        
        # 支出部分：該 Hub 的技術支援總成本
        - (c_hub[1] * d1)
    );

#Spoke 目標式
/*
minimize Spoke_obj:
         M_penalty * sum {j in 1..n}(Y_p[j] + Y_n[j]) 
       + delta * ( sum{i in 1..m} sum{j in 1..n-1} sum{k in j+1..n} (Z_p[i,j,k] + Z_n[i,j,k]))
       + sum{j in N} sum{i in M} ( (p[i]+ c_spoke[j]) * x[i,j] );

#spoke 選擇不同的hub有不同的p
*/
#hub 層級限制式

# 1. Budget constraint for each HUB Hospital
subject to budget_constraint_hub :
     c_hub[1] * d1 <= B_hub[1] + gamma* s[1]* sum{j in 1..n} x[1,j]; #c_hub[1], s[1]>>僅考慮Hub 1

# 2. ratio constraint for each HUB Hospital
subject to ratio_constraints_hub :
	p1 >= theta1 * c_hub[1]-theta2 * d1;

#KKT condition of Spoke(Spoke 層級限制式)

#FOC
s.t. Stationarity_x {i in 1..m, j in 1..n}:
     lam[i] + (mu[j]+ 1) * ( (if i==1 then p1 else p2) + c_spoke[j]) + pref_noise[i,j] + L[j] 
     - sum {k in j+1..n} R[i,j,k] + sum {k in 1..j-1} R[i,k,j] - xd[i,j] = 0; #納入對手價格

s.t. StationarityForY_p{j in 1..n}:
	M_penalty + L[j] - Yd_p[j] = 0;
    
s.t. StationarityForY_n{j in 1..n}:
	M_penalty - L[j] - Yd_n[j] = 0;
	
s.t. StationarityForZ_p{i in 1..m, j in 1..n-1, k in j+1..n}:
	delta + R[i,j,k] - Zd_p[i,j,k] = 0;
	
s.t. StationarityForZ_n{i in 1..m, j in 1..n-1, k in j+1..n}:
	delta - R[i,j,k] - Zd_n[i,j,k] = 0;


#Complementary Slackness >> Indicator Constraints
 # --- Capacity 互補 (lam * S_c = 0) ---
s.t. Ind_cap_1 {i in 1..m}: D_cap[i] == 0 ==> lam[i] == 0;
s.t. Ind_cap_2 {i in 1..m}: D_cap[i] == 1 ==> S_c[i] == 0;

# --- Budget 互補 (mu * S_b = 0) ---
s.t. Ind_budget_1 {j in 1..n}: D_budget[j] == 0 ==> mu[j] == 0;
s.t. Ind_budget_2 {j in 1..n}: D_budget[j] == 1 ==> S_b[j] == 0;

# --- x 變數互補 (x * xd = 0) ---
s.t. Ind_x_1 {i in 1..m, j in 1..n}: D_x[i,j] == 0 ==> x[i,j] == 0;
s.t. Ind_x_2 {i in 1..m, j in 1..n}: D_x[i,j] == 1 ==> xd[i,j] == 0;

# --- Y_p 變數互補 (Y_p * Yd_p = 0) ---
s.t. Ind_yp_1 {j in 1..n}: D_yp[j] == 0 ==> Y_p[j] == 0;
s.t. Ind_yp_2 {j in 1..n}: D_yp[j] == 1 ==> Yd_p[j] == 0;

# --- Z_p 變數互補 (Z_p * Zd_p = 0) ---
s.t. Ind_zp_1 {i in 1..m, j in 1..n-1, k in j+1..n}: D_zp[i,j,k] == 0 ==> Z_p[i,j,k] == 0;
s.t. Ind_zp_2 {i in 1..m, j in 1..n-1, k in j+1..n}: D_zp[i,j,k] == 1 ==> Zd_p[i,j,k] == 0;

# --- Z_n 變數互補 (Z_n * Zd_n = 0) ---
s.t. Ind_zn_1 {i in 1..m, j in 1..n-1, k in j+1..n}: D_zn[i,j,k] == 0 ==> Z_n[i,j,k] == 0;
s.t. Ind_zn_2 {i in 1..m, j in 1..n-1, k in j+1..n}: D_zn[i,j,k] == 1 ==> Zd_n[i,j,k] == 0;

#Constraints of Spoke
s.t. capacity_constraint_slack {i in 1..m}:
    sum {j in 1..n} x[i,j] + S_c[i] = (if i==1 then d1 else d2);
    #兩個hub的容量都納入考慮

s.t. cost_constraint_slack {j in 1..n}:
    sum {i in 1..m} ( (if i==1 then p1 else p2) + c_spoke[j]) * x[i,j] + S_b[j] = B_spoke[j];
    #兩個hub的價格都納入考量
    
s.t. y_var_constraint{j in 1..n}:
	Y_p[j] - Y_n[j] = T[j] - sum{i in 1..m}x[i,j];


s.t. z_var_constraint {i in 1..m, j in 1..n-1, k in j+1..n}:
    Z_p[i,j,k] - Z_n[i,j,k] = x[i,k] - x[i,j] ; 
