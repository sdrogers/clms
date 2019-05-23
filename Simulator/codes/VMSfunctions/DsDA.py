import glob
import pandas as pd


def get_schedule(n, schedule_dir):
    while True:
        files = sorted(glob.glob(schedule_dir + '\\*.csv'))
        if len(files) == n:
            last_file = files[-1]
            try:
                schedule = pd.read_csv(last_file)
                if schedule.shape[0] == 11951:
                    return last_file
            except:
                pass