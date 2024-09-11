import pandas as pd
import os
besler_data_path = r"C:\Users\hit\Desktop\liserSound\Data\BeslerData"
excel_file = r"C:\Users\hit\Desktop\liserSound\Data\BeslerData\תיק ניסוי מסודר.xlsx"
df = pd.read_excel(excel_file)
# Count occurrences of unique combinations of Tester and Word
unique_combinations = df[['Tester', 'Word']].drop_duplicates()
# Count occurrences of unique combinations of Tester and Word
combination_counts = df.groupby(['Tester', 'Word']).size().reset_index(name='Count')
countTotal=0
# Merge size count with unique combinations
unique_combinations2 = pd.merge(unique_combinations, combination_counts, on=['Tester', 'Word'], how='left')
countArray = unique_combinations2['Count']


for counting in countArray:
    for count_counting in range(counting):
        dfRow=df.iloc[count_counting+countTotal]
        tester = dfRow['Tester']
        date = dfRow['Date']
        word = str(dfRow['Word'])
        recording_number = dfRow['Recording Number']
        # Sort files by date
        tester_folder = os.path.join(besler_data_path, tester)
        word_folder = os.path.join(tester_folder, word)
        if count_counting==0:
            avi_files = [(filename, os.path.getctime(os.path.join(word_folder, filename))) for filename in os.listdir(word_folder) if filename.endswith(".avi")]
            avi_files.sort(key=lambda x: x[1])
            num_avi_files = len(avi_files) 

        if num_avi_files != counting:
            print(f"Error: Number of AVI files in folder ({tester}, {word}) does not match the recordings number")
            print(f"There are {num_avi_files} avi files and {counting} recordings")
            print("")
            continue
        # Find the file to rename
        file_count = avi_files[count_counting][0]
        old_file_path = os.path.join(word_folder, file_count)
        new_file_path = os.path.join(word_folder, f'{recording_number}.avi')
        if old_file_path != new_file_path:
            os.rename(old_file_path, new_file_path)
    countTotal += counting
