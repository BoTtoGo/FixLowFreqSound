from PyQt5.QtWidgets import QMainWindow, QApplication, QFileDialog, QLabel, QVBoxLayout
from PyQt5 import uic
import sys
import os
import librosa
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
import numpy as np
from scipy.io import wavfile
from scipy.signal import spectrogram
from scipy.signal import butter, filtfilt, lfilter
import torchaudio
import torch
from CleanUNet.network import CleanUNet

class MplCanvas(FigureCanvas):

    def __init__(self, parent=None):
        fig, self.ax = plt.subplots(3, 1, figsize=(8, 9))
        fig.subplots_adjust(hspace=1)
        super(MplCanvas, self).__init__(fig)


class UI(QMainWindow):
    def __init__(self):
        super(UI, self).__init__()

        # Load the ui file
        uic.loadUi("clean_sound_gui.ui", self)

        self.label_6 = self.findChild(QLabel, "label_6")
        self.label_8 = self.findChild(QLabel, "label_8")

        # Connect the buttons to functions
        self.pushButton.clicked.connect(self.browse_input)
        self.pushButton_4.clicked.connect(self.browse_output)
        self.pushButton_3.clicked.connect(lambda: self.denoise_audio(device='cpu'))

        # Initialize paths
        self.input_path = None
        self.output_path = None

        # Set up canvas for plotting
        self.input_canvas = MplCanvas(self)
        self.output_canvas = MplCanvas(self)
        
        input_layout = QVBoxLayout(self.frame)
        input_layout.addWidget(self.input_canvas)
        
        output_layout = QVBoxLayout(self.frame_3)
        output_layout.addWidget(self.output_canvas)

        # Show the app
        self.show()

    def browse_input(self):
        options = QFileDialog.Options()
        options |= QFileDialog.ReadOnly
        file_name, _ = QFileDialog.getOpenFileName(self, "Select Input WAV File", "", "WAV Files (*.wav)", options=options)
        if file_name:
            self.input_path = file_name
            self.textEdit.setText(file_name)

    def browse_output(self):
        options = QFileDialog.Options()
        options |= QFileDialog.ShowDirsOnly
        directory = QFileDialog.getExistingDirectory(self, "Select Output Directory", options=options)
        if directory:
            self.output_path = directory
            self.textEdit_3.setText(directory)


    # Function to load the trained model
    def load_model(self,model_path, device='cpu'):
        model = CleanUNet()
        model_state_dict = torch.load(model_path, map_location=torch.device(device), weights_only=True)
        model.load_state_dict(model_state_dict)
        model = model.to(device)
        model.eval()
        return model

    # Function to denoise a single audio file
    def denoise_audio(self, device='cpu'):
        output_audio_path = self.output_path
        noisy_audio_path = self.input_path
        model = self.load_model(r"cleanunet_model_conc.pth")

        noisy, sample_rate = torchaudio.load(noisy_audio_path)
        noisy = noisy.to(device)

        # Make sure the model and input are on the same device
        model = model.to(device)
        noisy = noisy.to(device)

        # Add batch dimension and pass through the model
        with torch.no_grad():
            denoised = model(noisy.unsqueeze(0))

        # Remove batch dimension and move back to CPU
        denoised = denoised.squeeze(0).cpu()
        output_audio_path = os.path.join(output_audio_path, "cleand.wav")
        # Save the denoised audio
        torchaudio.save(output_audio_path, denoised, sample_rate)

        self.plot_waveforms(noisy_audio_path, output_audio_path)


    def clean_sound(self):
        if not self.input_path or not self.output_path:
            print("Please select both input file and output path.")
            return

            # Read input file
        sample_rate, data = wavfile.read(self.input_path)

        # If the audio is stereo, take only one channel
        if data.ndim > 1:
            data = data[:, 0]

        # Design the band-pass filter
        lowcut = 1  # Low cutoff frequency in Hz
        highcut = 500  # High cutoff frequency in Hz
        nyquist = 0.5 * sample_rate
        low = lowcut / nyquist
        high = highcut / nyquist
        b, a = butter(3, [low, high], btype='band')

        # Apply the band-pass filter
        filtered_data = lfilter(b, a, data)

        # Perform the Fourier Transform on the filtered data
        n = len(filtered_data)
        frequency = np.fft.fftfreq(n, d=1 / sample_rate)
        magnitude = np.abs(np.fft.fft(filtered_data))

        # Save filtered data to output file
        output_file = f"{self.output_path}/cleaned_output.wav"
        wavfile.write(output_file, sample_rate, np.int16(filtered_data * 32768))

        # Plot input and output data
        self.plot_waveforms(self.input_path, output_file, sample_rate)

    def plot_waveforms(self, original_path, filtered_path):
        # Load audio files with librosa
        audio_original, sr1 = librosa.load(original_path, sr=None)
        audio_filtered, sr2 = librosa.load(filtered_path, sr=None)

        # Ensure the sample rates are the same for both audio files
        if sr1 != sr2:
            raise ValueError("Sample rates of original and filtered audio must be the same.")

        # Calculate the spectrogram for the original and filtered audio
        frequencies1, times1, Sxx1 = spectrogram(audio_original, fs=sr1, nperseg=1024, noverlap=512)
        Sxx1_db = 10 * np.log10(Sxx1 + np.finfo(float).eps)

        frequencies2, times2, Sxx2 = spectrogram(audio_filtered, fs=sr2, nperseg=1024, noverlap=512)
        Sxx2_db = 10 * np.log10(Sxx2 + np.finfo(float).eps)

        # Create time axis for waveform plots
        time_original = np.arange(audio_original.size) / sr1
        time_filtered = np.arange(audio_filtered.size) / sr2

        # Convert to dB scale for better visualization
        original_dB = 20 * np.log10(np.abs(audio_original) + np.finfo(float).eps)
        filtered_dB = 20 * np.log10(np.abs(audio_filtered) + np.finfo(float).eps)

        # Plot original waveform
        self.input_canvas.ax[0].plot(time_original, original_dB)
        self.input_canvas.ax[0].set_title('Input Waveform')
        self.input_canvas.ax[0].set_xlabel('Time [s]')
        self.input_canvas.ax[0].set_ylabel('Amplitude [dB]')

        # Create a dummy plot to attach an empty colorbar
        dummy_plot = self.input_canvas.ax[0].scatter([], [], c=[], cmap='Greys')
        self.input_canvas.figure.colorbar(dummy_plot, ax=self.input_canvas.ax[0])

        # Plot original spectrogram
        cax1 = self.input_canvas.ax[1].pcolormesh(times1, frequencies1, Sxx1_db, shading='gouraud', cmap='viridis')
        self.input_canvas.ax[1].set_title('Input Spectrogram')
        self.input_canvas.ax[1].set_xlabel('Time [s]')
        self.input_canvas.ax[1].set_ylabel('Frequency [Hz]')
        self.input_canvas.ax[1].set_ylim(0, 4000)  # Show only up to 5000 Hz
        # Add colorbar
        self.input_canvas.figure.colorbar(cax1, ax=self.input_canvas.ax[1], label='Intensity [dB]')

        # Plot original frequency spectrum
        frequency = np.fft.fftfreq(len(audio_original), d=1 / sr1)
        magnitude = 20 * np.log10(np.abs(np.fft.fft(audio_original)))
        self.input_canvas.ax[2].plot(frequency[:len(frequency) // 2], magnitude[:len(magnitude) // 2])
        self.input_canvas.ax[2].set_title('Input Frequency')
        self.input_canvas.ax[2].set_xlabel('Frequency [Hz]')
        self.input_canvas.ax[2].set_ylabel('Magnitude [dB]')
        self.input_canvas.ax[2].set_xlim(0, 4000)  # Show only up to 5000 Hz

        # Create a dummy plot to attach an empty colorbar
        dummy_plot_freq = self.input_canvas.ax[2].scatter([], [], c=[], cmap='Greys')
        self.input_canvas.figure.colorbar(dummy_plot_freq, ax=self.input_canvas.ax[2])

        # Plot filtered waveform
        self.output_canvas.ax[0].plot(time_filtered, filtered_dB)
        self.output_canvas.ax[0].set_title('Filtered Waveform')
        self.output_canvas.ax[0].set_xlabel('Time [s]')
        self.output_canvas.ax[0].set_ylabel('Amplitude [dB]')

        # Create a dummy plot to attach an empty colorbar
        dummy_plot_out = self.output_canvas.ax[0].scatter([], [], c=[], cmap='Greys')
        self.output_canvas.figure.colorbar(dummy_plot_out, ax=self.output_canvas.ax[0])

        # Plot filtered spectrogram
        cax2 = self.output_canvas.ax[1].pcolormesh(times2, frequencies2, Sxx2_db, shading='gouraud', cmap='viridis')
        self.output_canvas.ax[1].set_title('Filtered Spectrogram')
        self.output_canvas.ax[1].set_xlabel('Time [s]')
        self.output_canvas.ax[1].set_ylabel('Frequency [Hz]')
        self.output_canvas.ax[1].set_ylim(0, 4000)  # Show only up to 5000 Hz
        # Add colorbar
        self.output_canvas.figure.colorbar(cax2, ax=self.output_canvas.ax[1], label='Intensity [dB]')

        # Plot filtered frequency spectrum
        filtered_magnitude = 20 * np.log10(np.abs(np.fft.fft(audio_filtered)))
        self.output_canvas.ax[2].plot(frequency[:len(frequency) // 2],
                                      filtered_magnitude[:len(filtered_magnitude) // 2])
        self.output_canvas.ax[2].set_title('Filtered Frequency')
        self.output_canvas.ax[2].set_xlabel('Frequency [Hz]')
        self.output_canvas.ax[2].set_ylabel('Magnitude [dB]')
        self.output_canvas.ax[2].set_xlim(0, 4000)  # Show only up to 5000 Hz

        # Create a dummy plot to attach an empty colorbar
        dummy_plot_freq_out = self.output_canvas.ax[2].scatter([], [], c=[], cmap='Greys')
        self.output_canvas.figure.colorbar(dummy_plot_freq_out, ax=self.output_canvas.ax[2])

        self.input_canvas.draw()
        self.output_canvas.draw()

        self.label_6.setStyleSheet("background-color: #097099;")
        self.label_8.setStyleSheet("background-color: #097099;")


# Initialize The App
app = QApplication(sys.argv)
UIWindow = UI()
with open('SpyBot.qss', "r") as file:
    UIWindow.setStyleSheet(file.read())

app.exec_()
