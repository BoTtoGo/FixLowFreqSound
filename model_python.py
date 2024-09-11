import os
import tensorflow as tf
from keras_preprocessing import sequence
from scipy.io import wavfile
import numpy as np

def load_audio_files(noisy_dir, clean_dir):
    noisy_files = sorted(os.listdir(noisy_dir))
    clean_files = sorted(os.listdir(clean_dir))
    
    noisy_audio = []
    clean_audio = []
    
    for noisy_file, clean_file in zip(noisy_files, clean_files):
        sr_noisy, noisy_data = wavfile.read(os.path.join(noisy_dir, noisy_file))
        sr_clean, clean_data = wavfile.read(os.path.join(clean_dir, clean_file))
        
        if sr_noisy != sr_clean*6:
           raise ValueError("Sample rates do not match!")
        
        noisy_audio.append(noisy_data)
        clean_audio.append(clean_data)
    
    return np.array(noisy_audio), np.array(clean_audio)

noisy_dir=r"C:\Users\avimo\OneDrive\Desktop\NOISE"
clean_dir = r"C:\Users\avimo\OneDrive\Desktop\CLEAN"
noisy_audio, clean_audio = load_audio_files(noisy_dir, clean_dir)


def preprocess_audio(audio, max_length=None):
    if max_length:
        audio = sequence.pad_sequences(audio, maxlen=max_length, dtype='float32', padding='post', truncating='post')
    audio = np.expand_dims(audio, axis=-1)  # Adding channel dimension
    return audio

max_length = max(len(a) for a in noisy_audio)
noisy_audio = preprocess_audio(noisy_audio, max_length)
clean_audio = preprocess_audio(clean_audio, max_length)

from keras.models import Model
from keras.layers import Input, Conv1D, MaxPooling1D, UpSampling1D

def build_denoising_autoencoder(input_shape):
    input_audio = Input(shape=input_shape)

    # Encoder
    x = Conv1D(16, 3, activation='relu', padding='same')(input_audio)
    x = MaxPooling1D(2, padding='same')(x)
    x = Conv1D(8, 3, activation='relu', padding='same')(x)
    x = MaxPooling1D(2, padding='same')(x)
    x = Conv1D(8, 3, activation='relu', padding='same')(x)
    encoded = MaxPooling1D(2, padding='same')(x)

    # Decoder
    x = Conv1D(8, 3, activation='relu', padding='same')(encoded)
    x = UpSampling1D(2)(x)
    x = Conv1D(8, 3, activation='relu', padding='same')(x)
    x = UpSampling1D(2)(x)
    x = Conv1D(16, 3, activation='relu', padding='same')(x)
    x = UpSampling1D(2)(x)
    decoded = Conv1D(1, 3, activation='sigmoid', padding='same')(x)

    autoencoder = Model(input_audio, decoded)
    return autoencoder

input_shape = (max_length, 2)
autoencoder = build_denoising_autoencoder(input_shape)

autoencoder.compile(optimizer='adam', loss='mean_squared_error', metrics=['mae'])

history = autoencoder.fit(noisy_audio, clean_audio, epochs=50, batch_size=16, validation_split=0.2)
#---------

def denoise_audio(model, noisy_audio):
    denoised_audio = model.predict(np.expand_dims(noisy_audio, axis=0))
    return np.squeeze(denoised_audio)

test_file = 'path/to/test/noisy_audio.wav'
sr, noisy_test_audio = wavfile.read(test_file)
noisy_test_audio = preprocess_audio([noisy_test_audio], max_length)

denoised_audio = denoise_audio(autoencoder, noisy_test_audio[0])

# Save the denoised audio
wavfile.write('path/to/save/denoised_audio.wav', sr, denoised_audio)

#--------

import matplotlib.pyplot as plt

plt.plot(history.history['loss'], label='Training Loss')
plt.plot(history.history['val_loss'], label='Validation Loss')
plt.legend()
plt.show()


