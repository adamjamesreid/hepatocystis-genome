# This script is designed to take pseudobulk RNA-seq data e.g. MCA and generate mixtures at known percentages to test e.g. CIBERSORT for accuracy

import sys
import pandas as pd

# get file of pseudobulk
pb_file = sys.argv[1]
outstem = sys.argv[2]

pb = pd.read_csv(pb_file, sep='\t', index_col = 0)

mixtures = {
  'F20M20R20T20S20' :
    {'Female' : 0.2,
    'Male' : 0.2,
    'Schizont' : 0.2,
    'Trophozoite' : 0.2,
    'Ring' : 0.2},

  'F60M10R10T10S10' :
    {'Female' : 0.6,
    'Male' : 0.1,
    'Schizont' : 0.1,
    'Trophozoite' : 0.1,
    'Ring' : 0.1},

  'F10M60R10T10S10' :
    {'Female' : 0.1,
    'Male' : 0.6,
    'Schizont' : 0.1,
    'Trophozoite' : 0.1,
    'Ring' : 0.1},

  'F10M10R10T10S60' :
    {'Female' : 0.1,
    'Male' : 0.1,
    'Schizont' : 0.6,
    'Trophozoite' : 0.1,
    'Ring' : 0.1},

  'F25M25R25T25S0' :
    {'Female' : 0.25,
    'Male' : 0.25,
    'Schizont' : 0,
    'Trophozoite' : 0.25,
    'Ring' : 0.25}


}
# make mixtures
mix_df = pd.DataFrame(index=pb.index)
for m in mixtures:
  for s in mixtures[m]:
    if m not in mix_df:
      mix_df[m] = pb[s] * mixtures[m][s]
    else:
      mix_df[m] = mix_df[m] + (pb[s] * mixtures[m][s])
    
    #print(s, mix_df)

mix_df.to_csv(f'{outstem}.dat', sep='\t')
