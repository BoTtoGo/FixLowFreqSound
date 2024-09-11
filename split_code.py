import os
import numpy as np
import matplotlib.pyplot as plt
from pydub import AudioSegment
from pydub.silence import split_on_silence

def read_wav(audio_segment):
    """Convert AudioSegment to numpy array and get sample rate."""
    samples = np.array(audio_segment.get_array_of_samples())
    if audio_segment.channels == 2:
        samples = samples.reshape((-1, 2))
        samples = samples.mean(axis=1)
    return samples, audio_segment.frame_rate

def find_subarray_index(long_array, short_array):
    """Find the start index of short_array in long_array."""
    len_short = len(short_array)
    len_long = len(long_array)
    for i in range(len_long - len_short + 1):
        if np.array_equal(long_array[i:i+len_short], short_array):
            return i
    return -1

def calculate_dbfs(sample):
    """Calculate dBFS for a given sample."""
    return 20 * np.log10(np.abs(sample) / 32767.0 + 1e-6)  # add a small value to avoid log(0)

def calculate_silence_threshold(audio_segment, buffer_db=5):
    """Calculate the silence threshold dynamically."""
    return audio_segment.dBFS - buffer_db

microphone_dir = r"C:\Users\hit\Desktop\liserSound\Data\ModelData\afterSync_microphone"
laser_dir = r"C:\Users\hit\Desktop\liserSound\Data\ModelData\afterSync_Basler"
pos_option ='pos2' 
for mic_file in os.listdir(microphone_dir):
    if mic_file.endswith("_downsampled_aligned.wav"):
        base_name = os.path.splitext(mic_file)[0].replace("_downsampled_aligned", "")
        base_name = os.path.splitext(base_name)[0].replace("Recording_", "")
        laser_file = f"{base_name}_{pos_option}_upsampled_aligned.wav"
        
        mic_path = os.path.join(microphone_dir, mic_file)
        laser_path = os.path.join(laser_dir, laser_file)
        
        if not os.path.exists(laser_path):
            print(f"Matching file not found for {mic_file}")
            continue
        
        sound_file = AudioSegment.from_wav(mic_path)
        laser_sound_file = AudioSegment.from_wav(laser_path)

        # Get the sample rate
        sample_rate = sound_file.frame_rate
        laser_sample_rate = laser_sound_file.frame_rate

        # Convert the AudioSegment to raw audio data
        samples = np.array(sound_file.get_array_of_samples())

        # If stereo, take the mean of the two channels
        if sound_file.channels == 2:
            samples = samples.reshape((-1, 2))
            samples = samples.mean(axis=1)

        # Calculate dBFS for each sample
        dbfs = np.array([calculate_dbfs(sample) for sample in samples])

        # Create time array
        time = np.arange(len(samples)) / sample_rate

        # Plot the dBFS values over time
        plt.figure(figsize=(15, 5))
        plt.plot(time, dbfs)
        plt.title("Decibel Levels of the Audio File")
        plt.xlabel("Time (s)")
        plt.ylabel("Decibels (dBFS)")
        plt.show()

        # Calculate the silence threshold dynamically
        silence_thresh = calculate_silence_threshold(sound_file)
        print("Silence threshold (dBFS):", silence_thresh)

        # Split the audio file based on silence
        audio_chunks = split_on_silence(
            sound_file,
            min_silence_len=1000,
            silence_thresh=silence_thresh,
            keep_silence=1000
        )

        # Read the long wav file samples
        long_samples, long_frame_rate = read_wav(sound_file)

        # Process and export each chunk, and find its position in the long file
        for i, chunk in enumerate(audio_chunks):
            print(f"Chunk {i}:")
            mic_file = os.path.splitext(mic_file)[0].replace(".wav", "")
            out_file = os.path.join(microphone_dir, f"{mic_file}_{i}.wav")
            print("Exporting", out_file)
            chunk.export(out_file, format="wav")

            # Read the short wav file samples
            short_samples, short_frame_rate = read_wav(chunk)

            # Ensure both files have the same frame rate
            if long_frame_rate != short_frame_rate:
                raise ValueError("The frame rates of the two files do not match.")

            # Find the start index of the short file in the long file
            start_index = find_subarray_index(long_samples, short_samples)

            if start_index == -1:
                print("The chunk is not found in the long wav file.")
            else:
                # Calculate the start and end times
                start_time = start_index / long_frame_rate
                end_time = (start_index + len(short_samples)) / long_frame_rate

                print(f"  Start time: {start_time} seconds")
                print(f"  End time: {end_time} seconds")
                print(f"  Duration: {end_time - start_time} seconds")

                # Crop the corresponding segment from the laser recording
                laser_start_time_ms = start_time * 1000  # convert to milliseconds
                laser_end_time_ms = end_time * 1000  # convert to milliseconds

                laser_chunk = laser_sound_file[laser_start_time_ms:laser_end_time_ms]

                laser_file = os.path.splitext(laser_file)[0].replace(".wav", "")
                laser_out_file = os.path.join(laser_dir, f"{laser_file}_laser_{i}.wav")
                print("Exporting", laser_out_file)
                laser_chunk.export(laser_out_file, format="wav")