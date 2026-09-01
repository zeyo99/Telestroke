# Number of Hub Hospitals
param m;

# Number of spoke Hospitals
param n;

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

#Big M
param M_capacity;
param M_budget;
param BigM;

# Regularization parameter delta
param delta;

# Decision Variables
var x {i in 1..m, j in 1..n},integer >= 0; 
var Y_p{j in 1..n} >=0;
var Y_n{j in 1..n} >=0; 
var Z_p{i in 1..m, j in 1..n} >=0;
var Z_n{i in 1..m, j in 1..n} >=0;

#dual vars(對應primal之限制式)(亦已包含乘數非負限制)
var Lam{i in 1..m} >=0;
var mu{j in  1..n} >=0;
var V_l{j in 1..n} ;
var V_r{i in 1..m, j in 1..n} ;

var xd{i in 1..m, j in 1..n} >= 0;
var Yd_p{j in 1..n} >=0;
var Yd_n{j in 1..n} >=0; 
var Zd_p{i in 1..m, j in 1..n} >=0;
var Zd_n{i in 1..m, j in 1..n} >=0;

#binary vars(轉MILP)
var D_Lam{i in 1..m}, binary;
var D_mu{j in 1..n}, binary;
var D_Vl{j in 1..n}, binary;
var D_Vr{i in 1..m, j in 1..n}, binary;

# Objective Function: 
                              
minimize total_obj: #no sum fix
      ( sum {j in 1..n}(Y_p[j]+Y_n[j])
       + delta * sum {i in 1..m}sum{j in 1..n} (Z_p[i,j] + Z_n[i,j] ));

#FOC條件
s.t. Stationarity_x {i in 1..m, j in 1..n}:
    Lam[i]+ mu[j]* (p[i]+c_spoke[j]) + V_l[j] + (1-n)/n* V_r[i,j] - xd[i,j] = 0;

s.t. StationarityForY_p{j in 1..n}:
	1 + V_l[j] - Yd_p[j] = 0;
    
s.t. StationarityForY_n{j in 1..n}:
	1 - V_l[j] - Yd_n[j] = 0;
	
s.t. StationarutyForZ_p{i in 1..m, j in 1..n}:
	delta + V_r[i,j] - Zd_p[i,j] = 0;
	
s.t. StationarutyForZ_n{i in 1..m, j in 1..n}:
	delta - V_r[i,j] - Zd_n[i,j] = 0;
	
#Complementary Slackness

##capacity
s.t. Com_capacity{i in 1..m}:
	Lam[i] * (sum{j in 1..n}x[i,j] - d[i]) = 0;

##budget	
s.t. Com_budget{j in 1..n} :
	mu[j] * (sum{i in 1..m}(p[i] + c_spoke[j])* x[i,j] -B_spoke[j]) = 0;
	
##obj_left
s.t. Com_obj_left{j in 1..n}:
	V_l[j] * ((Y_p[j] - Y_n[j])-(T[j] - sum{i in 1..m}x[i,j])) = 0;
	
##obj_right
s.t. Com_obj_right{i in 1..m,j in 1..n}:
	V_r[i,j] * ((Z_p[i,j] - Z_n[i,j])-(x[i,j] - sum{h in 1..m, k in 1..n}x[h,k]/n)) = 0;

#Complementary Slackness-> Big M

#capacity
s.t. Com_capacity_1{i in 1..m}:
	sum{j in 1..n}x[i,j] - d[i] >= -M_capacity * D_Lam[i];

s.t. Com_capacity_2{i in 1..m}:
	Lam[i] <= M_capacity * (1- D_Lam[i]);
	
#budget
s.t. Com_budget_1{j in 1..n}:
	sum{i in 1..m}(p[i] + c_spoke[j])* x[i,j] -B_spoke[j] >= -M_budget * D_mu[j];
	
s.t. Com_budget_2{j in 1..n}:
	mu[j] <= M_budget * (1- D_mu[j]);

#obj_left
s.t. Com_obj_left_1{j in 1..n}:
	(Y_p[j] - Y_n[j])-(T[j] - sum{i in 1..m}x[i,j]) >= - BigM * D_Vl[j];

s.t. Com_obj_left_2{j in 1..n}:
	V_l[j] <= BigM * (1- D_Vl[j]);

#obj_right
s.t. Com_obj_right_1{i in 1..m,j in 1..n}:
	(Z_p[i,j] - Z_n[i,j])-(x[i,j] - sum{h in 1..m, k in 1..n}x[h,k]/n) >= - BigM * D_Vr[i,j];
	
s.t. Com_obj_right_2{i in 1..m,j in 1..n}:
	V_r[i,j] <= BigM * (1- D_Vr[i,j]);

# Constraints:

s.t. capacity_constraint {i in 1..m}:
    sum {j in 1..n} x[i,j] <= d[i];   
		
s.t. cost_constraint {j in 1..n}:
    sum {i in 1..m} (p[i] + c_spoke[j]) * x[i,j] <= B_spoke[j]; 

s.t. y_var_constraint{j in 1..n}:
	Y_p[j] - Y_n[j] = T[j] - sum{i in 1..m}x[i,j];
	
s.t. z_var_constraint{i in 1..m,j in 1..n}:
	Z_p[i,j] - Z_n[i,j] = x[i,j] - sum{h in 1..m, k in 1..n}x[h,k]/n;


