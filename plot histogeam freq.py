import numpy as np
from scipy.io import wavfile
import matplotlib.pyplot as plt
import wave
import struct

# Step 1: Read the WAV file
wav_file_path = r"C:\Users\avimo\OneDrive\Desktop\ModelData\originals files\Recording (682).wav"


# Read the WAV file
with wave.open(wav_file_path, 'r') as wav_file:
    # Extract raw audio from WAV file
    sample_rate = wav_file.getframerate()
    num_frames = wav_file.getnframes()
    raw_data = wav_file.readframes(num_frames)
    num_channels = wav_file.getnchannels()

# Convert raw data to numpy array
fmt = f"{num_frames * num_channels}h"
data = struct.unpack(fmt, raw_data)
data = np.array(data)

# If stereo, take one channel (for simplicity)
if num_channels > 1:
    data = data[::num_channels]

# Perform Fourier Transform
frequencies = np.fft.rfftfreq(len(data), 1 / sample_rate)
magnitude = np.abs(np.fft.rfft(data))

# Define frequency bands (similar to your provided histogram)
freq_bands = [125, 250, 500, 1000, 2000, 4000, 8000]
intelligibility = []

# Calculate the percentage of total magnitude for each band
for i in range(len(freq_bands) - 1):
    band = (frequencies >= freq_bands[i]) & (frequencies < freq_bands[i + 1])
    band_magnitude = np.sum(magnitude[band])
    total_magnitude = np.sum(magnitude)
    intelligibility.append((band_magnitude / total_magnitude) * 100)

# Append the last frequency band
band = frequencies >= freq_bands[-1]
band_magnitude = np.sum(magnitude[band])
intelligibility.append((band_magnitude / total_magnitude) * 100)

# Plot the histogram
plt.bar([str(band) for band in freq_bands], intelligibility)
plt.xlabel('Frequency (Hz)')
plt.ylabel('Intelligibility (%)')
plt.title('Frequency Intelligibility Histogram')
plt.grid(True)
plt.show()