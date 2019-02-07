import pandas as pd

SEQ_HEAVY_IC_CUTOFF = 1.0

df = pd.read_csv('../results/d0/publish/denriched_same.4.l1.r1.csv')

print(df.head())
print(sum(df['tf.x'] != df['tf.y']))


df = df[['kmer', 'label', 'tf.x', 'primer.x', 'family.x', 'en_th', 'shape']].rename(index=str, columns={'kmer':'shapemer', 'tf.x':'tf', 'primer.x':'primer', 'family.x':'family'})

print(df.head())

mics = pd.read_csv('../results/d0/publish/median_ics.4.l1.r1.csv')
mics = mics[mics['medianic'] > SEQ_HEAVY_IC_CUTOFF]

print(mics)

df = df.merge(mics)

fsbs = lambda x: sum(x == 'both')
fsbs.__name__ = 'FSBSIgnore'


fnbs = lambda x: sum(x == 'bg')
fnbs.__name__ = 'FNBSIgnore'

df = df.groupby(['en_th', 'shape', 'family', 'tf', 'primer']).agg({'label': [fsbs, fnbs]})
df.columns = df.columns.droplevel(0)
df = df.reset_index()

df.to_csv('../results/d0/publish/mic_correction_same.4.l1.r1.csv', index=False)



