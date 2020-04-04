import pandas as pd
import numpy as np

# df = pd.read_csv('old_IVT/test_campbell_test_campbell_030_fixed.tsv',delimiter='\t')
df = pd.read_csv('old_IVT_MovingMedian/test_campbell_test_campbell_030_fixed.tsv',delimiter='\t')

df.index = pd.to_datetime(df['EyeTrackerTimestamp'], unit='us')

# print(df)

def wrapperfun(x):
    # print('fewwef')
    valArray = np.isnan(x)
    if any(valArray):
        return 'NaNs: ' + str(np.sum(valArray)) + '/' + str(len(valArray))
    else:
        return np.mean(x)

test = df.groupby(pd.Grouper(freq='500ms'))['FixationPointX (MCSpx)', 'FixationPointY (MCSpx)'].aggregate([wrapperfun])

test.to_csv('temp2.csv')
print(test)