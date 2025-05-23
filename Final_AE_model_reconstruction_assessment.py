import numpy as np
from keras.models import Model
from keras.layers import Input, Dense, BatchNormalization, Dropout
from keras import regularizers
from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pickle
import os

import tensorflow as tf
from tensorflow.keras.models import Model
from tensorflow.keras.layers import Input, Dense, Dropout, BatchNormalization
from tensorflow.keras.callbacks import EarlyStopping, ReduceLROnPlateau
from tensorflow.keras import regularizers
from tensorflow.keras.models import load_model

from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error
from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score
from sklearn.metrics import mean_absolute_percentage_error
from sklearn.decomposition import PCA
#from sklearn.preprocessing import StandardScaler   #No need for further standiration, data is already normalized and log transformed  (No??)
from sklearn.preprocessing import MinMaxScaler
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score

Sigs = r'C:\Users\lyaziddy\Desktop\XAIM_PAVIA\AE\all_log1p_sigs_dea.tsv'      #CHANGEE !!!
df1 = pd.read_csv(Sigs, sep="\t")
df1.rename(columns={'Unnamed: 0': 'Hugo_Symbol'}, inplace=True)
df1_transposed = df1.set_index('Hugo_Symbol').transpose()
#df1_transposed.reset_index(inplace=True)
print(df1_transposed)
print(df1_transposed.index)


Meta = r'C:\Users\lyaziddy\Desktop\XAIM_PAVIA\AE\Metadata_IDs_PatientVSNormal_final.tsv'
df2 = pd.read_csv(Meta, sep="\t")



min_values = df1_transposed.min()    #minimum values
max_values = df1_transposed.max()    #maximum values
range_values = max_values - min_values    # Calculate the range
mean_values = np.mean(df1_transposed, axis=0)
std_values = np.std(df1_transposed, axis=0)
coefficient_of_variation = std_values / mean_values    # Calculate the coefficient of variation for each feature of the original (non-normalized) data
num_data_2 = pd.DataFrame({
    'Min': min_values,
    'Max': max_values,
    'Range': range_values,
    'CV': coefficient_of_variation
})
print("Scale of each feature (Min, Max, Range):\n", num_data_2)



#Remove irrelevant features (based on the CV)     #to help reduce the dim even before applying any dim red technique
#CV Should be computed on the original data, as it does not have a meaningful interpretation after standardization.
#step: Filter Before Standardization: Use the CV to filter your data before applying standardization.

filtered_num_data_2 = num_data_2[num_data_2['CV'] > 0.1]
print(filtered_num_data_2)       #STOPPED HERE
features_to_keep = filtered_num_data_2.index
filtered_df1_transposed = df1_transposed[features_to_keep.intersection(df1_transposed.columns)]
print(filtered_df1_transposed)      #[1187 rows x 4507 columns]       Reduced dim from 5169 to 4507 (with CV > 0.1)


# Normalize the data to range [0, 1]                 #MMinMaxScaler might NOT be the best option, When values are scaled to 0 to 1, very small values can become significant when calculating % errors like MAPE. This can lead to extremely high MAPE values.
#scaler = MinMaxScaler()
#num_data = df1_transposed.select_dtypes(include=[np.number])   # Separate the numerical columns  #No need, already set the patients IDs as index
scaler = StandardScaler()         #data will be standardized to have a mean of 0 and a standard deviation of 1
num_scaled = scaler.fit_transform(filtered_df1_transposed)   # Fit and transform the numerical data
scaled_df = pd.DataFrame(num_scaled, columns=filtered_df1_transposed.columns, index=filtered_df1_transposed.index)    # Convert the scaled data back to a DataFrame
print(scaled_df)


#Before running the PCA: split the scaled data into training(80%), test(10%) and validation(10%) based on class labels (Normal vs Patients)

scaled_df.reset_index(inplace=True)
#print(scaled_df)
merged_df = pd.merge(scaled_df, df2, left_on='index', right_on='Sample_ID', how='inner')
merged_df = merged_df.drop(columns=['Sample_ID', 'Patient_ID'])
merged_df = merged_df.set_index('index')
#print(merged_df)
#sklearn.model_selection.train_test_split(*arrays, test_size=None, train_size=None, random_state=None, shuffle=True, stratify=None)
scaled_train, scaled_X = train_test_split(merged_df, test_size=0.2, stratify=merged_df['Sample_type'])           #split the data into training (949), test (119) and validation (119)
scaled_test, scaled_val = train_test_split(scaled_X, test_size=0.5, stratify=scaled_X['Sample_type'])
#print(scaled_train)
#print(scaled_test)
#print(scaled_val)


# Normalize the data
#data = df1.values
#data_norm = (data - np.mean(data)) / np.std(data)

# Split data into train and test sets (or use all data for demonstration)
# Split data into training and testing sets
#x_train, x_test = train_test_split(data_scaled, test_size=0.2, random_state=42)       #If more data is available, we could have split also to validation set
# Further split the training data into training and validation sets
#x_train, x_val = train_test_split(x_train, test_size=0.2, random_state=42)    #PUT A STRATIFICATION BASED ON the sample type! DONE!

scaled_train = scaled_train.drop(columns=['Sample_type'])
scaled_test = scaled_test.drop(columns=['Sample_type'])
scaled_val = scaled_val.drop(columns=['Sample_type'])

# Convert data to float32 explicitly
x_train = np.array(scaled_train, dtype=np.float32)
x_test = np.array(scaled_test, dtype=np.float32)
x_val = np.array(scaled_val, dtype=np.float32)
#x_train = data_norm


# Define the custom R-squared function (from Alessandro)
def custom_r2(X, X_reconstr):
    X = np.asarray(X)
    X_reconstr = np.asarray(X_reconstr)
    return 1 - sum((X - X_reconstr) ** 2).sum() / sum((X - X.mean()) ** 2).sum()


# Define the performance evaluation function
def eval_perf(X, X_reconstr):
    # Mean Squared Error
    mse = mean_squared_error(X, X_reconstr)
    # Mean Absolute Error
    mae = mean_absolute_error(X, X_reconstr)
    # R-squared
    r2_sklearn = r2_score(X, X_reconstr, multioutput='uniform_average')
    # Custom R-squared
    r2_manual = custom_r2(X, X_reconstr)
    # Mean Absolute Percentage Error (MAPE)
    mape = np.mean(np.abs((X - X_reconstr) / X)) * 100
    return {'MSE': mse, 'MAE': mae, 'R2_sklearn': r2_sklearn, 'r2_custom': r2_manual, 'MAPE': mape}


# Define the autoencoder model
def AE_model(input_dim, latent_dim):
    input_layer = Input(shape=(input_dim,))
    # Encoder
    encoded = Dense(4000, activation='relu', kernel_regularizer=regularizers.l2(0.001))(input_layer)
    encoded = BatchNormalization()(encoded)
    encoded = Dropout(0.1)(encoded)

    encoded = Dense(3000, activation='relu', kernel_regularizer=regularizers.l2(0.01))(encoded)
    encoded = BatchNormalization()(encoded)
    encoded = Dropout(0.1)(encoded)

    encoded = Dense(2000, activation='relu', kernel_regularizer=regularizers.l2(0.001))(encoded)
    encoded = BatchNormalization()(encoded)
    encoded = Dropout(0.1)(encoded)

    encoded = Dense(1000, activation='relu', kernel_regularizer=regularizers.l2(0.001))(encoded)
    encoded = BatchNormalization()(encoded)
    encoded = Dropout(0.1)(encoded)

    # Latent Space
    encoded_output = Dense(latent_dim, activation='relu')(encoded)

    # Decoder
    decoded = Dense(1000, activation='relu', kernel_regularizer=regularizers.l2(0.001))(encoded_output)
    decoded = BatchNormalization()(decoded)
    decoded = Dropout(0.1)(decoded)

    decoded = Dense(2000, activation='relu', kernel_regularizer=regularizers.l2(0.001))(decoded)
    decoded = BatchNormalization()(decoded)
    decoded = Dropout(0.1)(decoded)

    decoded = Dense(3000, activation='relu', kernel_regularizer=regularizers.l2(0.001))(decoded)
    decoded = BatchNormalization()(decoded)
    decoded = Dropout(0.1)(decoded)

    decoded = Dense(4000, activation='relu', kernel_regularizer=regularizers.l2(0.001))(decoded)
    decoded = BatchNormalization()(decoded)
    decoded = Dropout(0.1)(decoded)

    decoded_output = Dense(input_dim, activation='linear')(decoded)

    encoder = Model(input_layer, encoded_output)
    autoencoder = Model(input_layer, decoded_output)

    return autoencoder, encoder


# Specify the input dimension
input_dim = 4507  # Adjust this according to your input shape

# Specify latent dimensions for the models
latent_dims = [100, 200, 449]

# Evaluate reconstruction metrics for each latent dimension
for latent_dim in latent_dims:
    # Build the autoencoder model
    autoencoder, encoder = AE_model(input_dim, latent_dim)

    # Load the saved weights
    weights_path = f'C:/Users/lyaziddy/Desktop/XAIM_PAVIA/AE/model_weights_{latent_dim}.weights.h5'
    autoencoder.load_weights(weights_path)

    # Prepare the training data
    x_train = np.array(scaled_train, dtype=np.float32)

    # Make predictions
    x_train_pred = autoencoder.predict(x_train)

    # Evaluate performance
    perf_metrics = eval_perf(x_train, x_train_pred)

    # Output the metrics
    print(f'Reconstruction metrics for latent dimension {latent_dim}: {perf_metrics}')
