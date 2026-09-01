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

# Decision Variables
var x {i in M, j in N} >= 0; 
var Y_p{j in N} >= 0;
var Y_n{j in N} >= 0; 
var Z_p{i in 1..m, j in 1..n-1, k in j+1..n} >= 0;
var Z_n{i in 1..m, j in 1..n-1, k in j+1..n} >= 0;
# Objective Function: 
                              
minimize total_obj:
       M_penalty * sum {j in 1..n}(Y_p[j] + Y_n[j]) 
       + delta * ( sum{i in 1..m} sum{j in 1..n-1} sum{k in j+1..n} (Z_p[i,j,k] + Z_n[i,j,k]))
       + sum{j in N} sum{i in M} ( (p[i]+ c_spoke[j]) * x[i,j] ); 

# Constraints:

s.t. capacity_constraint {i in 1..m}:
    sum {j in 1..n} x[i,j] <= d[i];   #lam

s.t. cost_constraint {j in 1..n}:
    sum {i in 1..m} ( p[i] + c_spoke[j]) * x[i,j] <= B_spoke[j]; #mu

s.t. y_var_constraint{j in 1..n}:
	Y_p[j] - Y_n[j] = T[j] - sum{i in 1..m}x[i,j]; #L

s.t. z_var_constraint{i in 1..m, j in 1..n-1, k in j+1..n}:
	Z_p[i,j,k] - Z_n[i,j,k] = x[i,k] - x[i,j]; #R






