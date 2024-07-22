import pandas as pd
import matplotlib.pyplot as plt

# Read the CSV files
prod_costs_df = pd.read_csv('output/Prod_Costs.csv', header=None)
prices_0_df = pd.read_csv('output/prices_0.csv', header=None)
quants_0_df = pd.read_csv('output/quants_0.csv', header=None)

# Ensure the lengths match
assert len(prod_costs_df) == len(prices_0_df) == len(quants_0_df), "The lengths of the CSV files do not match."

# Calculate the vertices of the parabolas and scale them
vertices = []

for index, row in prod_costs_df.iterrows():
    c1, c2, c3 = row[1], row[2], row[3] 
    # Vertex formula: x = -b / (2a) for ax^2 + bx + c
    # Here, a = 3*c3, b = 2*c2
    x_vertex = -c2 / (3 * c3)
    y_vertex = c1 - (c2**2) / (3 * c3)

    # Scale the vertices
    x_vertex_scaled = x_vertex / quants_0_df.iloc[index, 0]
    y_vertex_scaled = y_vertex / prices_0_df.iloc[index, 0]

    vertices.append((x_vertex_scaled, y_vertex_scaled))

# Convert vertices to DataFrame for easier plotting
vertices_df = pd.DataFrame(vertices, columns=['x', 'y'])

# Plot the vertices
plt.scatter(vertices_df['x'], vertices_df['y'])
plt.title('Scaled Vertices of Parabolas')
plt.xlabel('x (percentage of quants_0)')
plt.ylabel('y (percentage of prices_0)')
plt.grid(True)

# Plot the line with x- and y-intercepts both 1
x_vals = [0, 1]
y_vals = [1, 0]
plt.plot(x_vals, y_vals, 'r--')

# Save the plot as an image file
plt.savefig('output/scaled_vertices_plot.png')
