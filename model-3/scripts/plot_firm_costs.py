import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys

def main(k):
    # Read values from CSV files
    q0 = pd.read_csv('output/quants_0.csv', header=None).iloc[k, 0]
    p0 = pd.read_csv('output/prices_0.csv', header=None).iloc[k, 0]
    ela_ii = pd.read_csv('output/Elas.csv', header=None).iloc[k, k] 
    prod_costs = pd.read_csv('output/Prod_Costs.csv', header=None).iloc[k, :]
    c0, c1, c2, c3 = prod_costs[0], prod_costs[1], prod_costs[2], prod_costs[3]

    q = np.linspace(1, q0, 100)
    q2 = np.linspace(q0 / 10, q0, 90)
    demand = p0 - (p0 / q0) * q
    marginal_cost = c1 + 2 * c2 * q + 3 * c3 * q**2
    avg_total_cost = c0/q2 + c1 + c2*q2 + c3*q2**2
    avg_var_cost = c1 + c2*q + c3*q**2
    marginal_revenue = p0 - (p0/q0 + 1/ela_ii) * q   # does not take into account cross-price elasticities 

    # Plotting
    plt.figure(figsize=(10, 6))
    plt.plot(q, demand, label='demand', color='navy')
    plt.plot(q, marginal_revenue, label='marginal revenue', color='blue') 
    plt.plot(q, marginal_cost, label='marginal cost', color='red')
    plt.plot(q2, avg_total_cost, label='ATC', color='orange')
    plt.plot(q, avg_var_cost, label='AVC', color='gold')

    # Adding labels and title
    plt.xlabel('quantity')
    plt.ylabel('USD')
    plt.xlim((0, q0)) 
    plt.ylim((0, p0 * 1.2))
    plt.title(f'Costs & Revenue for firm {k}')
    plt.axhline(0, color='black', linewidth=0.5)
    plt.axvline(0, color='black', linewidth=0.5)
    plt.grid(color='gray', linestyle='--', linewidth=0.5)
    plt.legend()
    plt.savefig('output/firm_costs.png')


if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("Usage: python script.py <k>")
        sys.exit(1)
    
    k = int(sys.argv[1])
    main(k)
