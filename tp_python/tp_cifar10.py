import numpy as np
import tensorflow as tf
from keras.datasets import cifar10

#mnist = tf.keras.datasets.mnist

(x_train, y_train),(x_test, y_test) = cifar10.load_data()
print('x_train shape : ' , x_train.shape)
x_train, x_test = x_train /255.0, x_test/255.0

model = tf.keras.models.Sequential([
    tf.keras.layers.Flatten(input_shape = [32,32,3]),
    tf.keras.layers.Dense(512,activation = 'relu'),
    tf.keras.layers.Dense(1024,activation = 'relu'),
    tf.keras.layers.Dense(512,activation = 'relu'),
    tf.keras.layers.Dense(10, activation = 'softmax')
])

model.compile(
    optimizer= 'adam',
    loss = 'sparse_categorical_crossentropy',
    metrics = ['accuracy']
)

model.summary()

model.fit(x_train, y_train, epochs = 10)
model.evaluate(x_test,y_test,verbose = 2)
