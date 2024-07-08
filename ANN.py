import pandas as pd
from sklearn.preprocessing import MinMaxScaler
from sklearn.model_selection import train_test_split
import numpy as np
import tensorflow as tf
from tensorflow.keras.models import Model
from tensorflow.keras.layers import Input, Dense
import matplotlib.pyplot as plt

# Load the SIDRstats table
file_path = './data/ogSIDRstats.tsv'
sidr_stats_df = pd.read_csv(file_path, sep='\t')

# Select numerical columns
numerical_cols = [
    'Avg_fold', 'Length', 'Ref_GC', 'Covered_bases', 'Covered_percent',
    'Read_GC', 'Plus_reads', 'Minus_reads', 'RNA_Covered_bases', 'RNA_Avg_fold',
    'RNA_Covered_percent', 'RNA_Read_GC'
]

# Initialize the scaler
scaler = MinMaxScaler()

# Fit and transform the numerical columns
sidr_stats_df[numerical_cols] = scaler.fit_transform(sidr_stats_df[numerical_cols])

# Filter the data for 'Oscheius dolichura'
true_target_df = sidr_stats_df[sidr_stats_df['Origin'] == 'Oscheius dolichura']

# Use only 1/3 of the data for training
train_df = true_target_df.sample(frac=1/3, random_state=42)

# Extract numerical data for training
X = train_df[numerical_cols].values

# Function to create and train an autoencoder
def create_autoencoder(input_dim, encoding_dim):
    input_layer = Input(shape=(input_dim,))
    encoded = Dense(encoding_dim, activation='relu')(input_layer)
    decoded = Dense(input_dim, activation='sigmoid')(encoded)
    autoencoder = Model(inputs=input_layer, outputs=decoded)
    autoencoder.compile(optimizer='adam', loss='mean_squared_error')
    return autoencoder


# Define a range of encoding dimensions to test
encoding_dims = [2, 5, 10, 20, 50]
errors = []

# Split the data into training and validation sets
X_train, X_val = train_test_split(X, test_size=0.2, random_state=42)

for encoding_dim in encoding_dims:
    print(f'Training autoencoder with encoding dimension: {encoding_dim}')
    autoencoder = create_autoencoder(input_dim=X.shape[1], encoding_dim=encoding_dim)
    autoencoder.fit(X_train, X_train, epochs=50, batch_size=32, shuffle=True, validation_data=(X_val, X_val), verbose=0)

    # Evaluate the model on the validation set
    encoded_data = autoencoder.predict(X_val)
    reconstruction_error = np.mean((X_val - encoded_data) ** 2, axis=1)
    mean_error = np.mean(reconstruction_error)
    errors.append(mean_error)
    print(f'Mean Reconstruction Error for encoding_dim {encoding_dim}: {mean_error}')

# Find the encoding dimension with the lowest error
best_encoding_dim = encoding_dims[np.argmin(errors)]
print(f'Best encoding dimension: {best_encoding_dim}')

# Plot the errors
plt.plot(encoding_dims, errors, marker='o')
plt.xlabel('Encoding Dimension')
plt.ylabel('Mean Reconstruction Error')
plt.title('Reconstruction Error vs. Encoding Dimension')
plt.show()

print(f'Errors for each encoding dimension: {errors}')
print(f'Best encoding dimension: {best_encoding_dim}')

######################################################
# Running the Autoencoder on the Full Dataset #######
#####################################################

# Create and train the autoencoder
autoencoder = create_autoencoder(input_dim=X_train.shape[1], encoding_dim=best_encoding_dim)
autoencoder.fit(X_train, X_train, epochs=50, batch_size=32, shuffle=True, validation_split=0.2, verbose=0)

# Use the trained autoencoder to reconstruct the full dataset
X_full = sidr_stats_df[numerical_cols].values
encoded_data = autoencoder.predict(X_full)
reconstruction_error = np.mean((X_full - encoded_data) ** 2, axis=1)

# Define the 5% reconstruction error threshold
threshold = 0.05

# Categorize the contigs based on the reconstruction error
sidr_stats_df['Reconstruction_Error'] = reconstruction_error
sidr_stats_df['Category'] = ['Normal' if error < threshold else 'Anomaly' for error in reconstruction_error]

# Separate the normal and anomaly contigs
normal_contigs = sidr_stats_df[sidr_stats_df['Category'] == 'Normal']
anomaly_contigs = sidr_stats_df[sidr_stats_df['Category'] == 'Anomaly']

# Save the results to CSV files
normal_contigs.to_csv('normal_contigs.csv', index=False)
anomaly_contigs.to_csv('anomaly_contigs.csv', index=False)

print(f'Normal contigs: {len(normal_contigs)}')
print(f'Anomaly contigs: {len(anomaly_contigs)}')