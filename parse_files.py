import pandas as pd
import numpy as np
import os,sys

path_prefix = sys.argv[1]
for filename in os.listdir(path_prefix):
    if filename.endswith(".tsv"):
        df = pd.read_csv(path_prefix+'/'+filename, delimiter='\t')
        df = df.iloc[:, :-1] # Due to trailing tab

        # where true, zero out values
        df['remove_info'] = False

        # Query info
        df['remove_info'] += ~(df['GazeEventType'].str.match('^Fixation'))
        df['remove_info'] += ~(df['ValidityLeft'] < 3)
        df['remove_info'] += ~(df['ValidityRight'] < 3)

        # Where true, NaN these: 'FixationPointX (MCSpx)', 'FixationPointY (MCSpx)'
        mask = (df['remove_info'] == True)
        df.loc[mask, 'FixationPointX (MCSpx)'] = np.NaN
        df.loc[mask, 'FixationPointY (MCSpx)'] = np.NaN

        df.drop('remove_info', axis=1, inplace=True)

        new_filename = filename.split('.')[0]
        new_filename += '_fixed.tsv'

        df.to_csv(new_filename,sep='\t',index=False) # ,na_rep=''