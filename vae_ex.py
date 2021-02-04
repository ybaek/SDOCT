"""Demo code to practice
   
   1. Importing R objects into Python
   2. Training TF-neural net structures
"""
import numpy as np
import pickle
import rpy2.robjects as robjects
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import backend as K
from tensorflow.keras import layers

# Importing R object + preprocessing
Z_DIM = 288
Y_DIM = 64
readRDS = robjects.r['readRDS']
rdso    = readRDS('~/projects/SDOCT/data/fit_data.Rds')
z = np.array(rdso[0])
y = np.array(rdso[1])
## Missing rows will be omitted
nan_entries = np.isnan(y)
nan_rows    = np.unique(np.where(nan_entries)[0])
z, y        = np.delete(z, nan_rows, axis = 0), np.delete(y, nan_rows, axis = 0)
## Linear interpolation map to (-1, 1)
z_s = np.interp(z, (z.min(), z.max()), (-1., +1.))
y_s = np.interp(y, (y.min(), y.max()), (-1., +1.))
## Dataset object
z_s = tf.data.Dataset.from_tensor_slices(z_s)
y_s = tf.data.Dataset.from_tensor_slices(y_s)
data_s = tf.data.Dataset.zip( (z_s, y_s) )
data_s = data_s.shuffle(1000)
## Train-test split
data_train = data_s.skip(200)
data_test  = data_s.take(200)

def _unzip_data(zipdata):
    X = np.reshape(np.concatenate([x for x,y in zipdata], axis = 0), (-1, Z_DIM))
    Y = np.reshape(np.concatenate([y for x,y in zipdata], axis = 0), (-1, Y_DIM))
    return (X, Y)

data_train = _unzip_data(data_train)
data_test  = _unzip_data(data_test)
z_train, y_train = data_train[0], data_train[1]
z_test, y_test   = data_test[0], data_test[1]
y_train, y_test = np.reshape(y_train, (-1, 8, 8)), np.reshape(y_test, (-1, 8, 8))
with open('testset.pkl', 'wb') as f:
    pickle.dump((z_test, y_test), f, protocol = pickle.HIGHEST_PROTOCOL)

"""
The following code for training the VAE is not yet incorporating the actual data
Code to construct a custom neural net architecture
"""
# Global constants
IN_DIM = (288, 1)
OUT_DIM = (8, 8, 1)
IN_LD = 10
OUT_LD = 8
# Model will consist of an "encoder" and a "decoder"
inputs    = keras.Input(shape = IN_DIM)
CL_1      = layers.Conv1D(1, 2, 2, activation = 'relu')(inputs)
CL_2      = layers.Conv1D(1, 2, 2, activation = 'relu')(CL_1)
RL_1      = layers.Flatten()(CL_2)
latent_x  = layers.Dense(IN_LD, activation = 'relu')(RL_1)
# Fully connected layer (Linear model: no activation!)
# Ridge penalty corresponds to a simple Gaussian prior
latent_y    = layers.Dense(OUT_LD,
                           kernel_regularizer = keras.regularizers.l2(1.0),
                           bias_regularizer = keras.regularizers.l2(1.0))(latent_x)
DL_1        = layers.Dense(4 * 4 * 1, activation = 'relu')(latent_y)
RL_2        = layers.Reshape((4, 4, 1))(DL_1)
outputs        = layers.Conv2DTranspose(1, 3, (2, 2),
                                     activation = 'sigmoid', padding = 'same')(RL_2)
model       = keras.Model(inputs, outputs, name = 'model')
model.summary()
# Compile and run
EPOCHS = 50
BATCH_SIZE = 5
model.compile(optimizer = 'adam',
              loss = keras.losses.mean_squared_error)
history = model.fit(z_train, y_train, 
                    batch_size = BATCH_SIZE, epochs = EPOCHS, 
                    validation_data = (z_test, y_test))
# Ridge loss print-out
yhat_test = np.array(model.predict(z_test))
model.save('model.h5')
with open('yhat.pkl', 'wb') as f:
    pickle.dump(yhat_test, f, protocol = pickle.HIGHEST_PROTOCOL)
