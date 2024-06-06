import matplotlib.pyplot as plt
import pandas as pd

# Sample data (replace this with your actual data)
data = {
    'Length': [100, 200, 150, 152, 250],
    'Coverage': [50, 60, 70, 71, 90],
    'GC_Content': [40, 45, 80, 55, 60],
    'Origin': ['A', 'B', 'C', 'D', 'E']
}

df = pd.DataFrame(data)

# Map 'Origin' to colors
color_dict = {'A': 'red', 'B': 'blue', 'C': 'green', 'D': 'orange', 'E': 'purple'}

# Plot
plt.scatter(df['Length'], df['Coverage'], s=df['GC_Content']*10, c=df['Origin'].map(color_dict), alpha=0.5)

# Add labels and title
plt.xlabel('Length')
plt.ylabel('Coverage')
plt.title('Bubble Plot')

# Show plot
plt.show()
