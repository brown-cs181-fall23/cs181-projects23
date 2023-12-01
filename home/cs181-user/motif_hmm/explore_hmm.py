import hmmlearn as hmm
import numpy as np
from preprocess import preprocess_motif
from hmm_traintest import train_hmm, test_hmms
import pickle


## DO NOT MODIFY
motifs = ["chinmo", "kenbarbie", "lameduck", "sugarbabe"]
bed_files = ["data/" + motif + "_motifs.fa" for motif in motifs]

preprocess_motif(motifs, bed_files)

train_data = "train_motif_dict.pkl"
test_data = "test_motif_dict.pkl"
mystery_bag = "data/mystery.pkl"
## DO NOT MODIFY**


# Train HMM
trained_models = {}  # key: motif, value: trained HMM
for motif in motifs:
    # TODO (Task 2d): fill in the training step.
    pass

# Test HMM
results = test_hmms(trained_models, test_data, True)

print("Results for TEST data.")
for motif, result in zip(motifs, results):
    print(f"{motif}: {result}")

print()

# UNCOMMENT IN FINAL STEP (exploring mystery bag data, Task 4)
# mystery_results = test_hmms(trained_models, mystery_bag, False)

# print("Results for mystery data.")
# for motif, result in zip(motifs, mystery_results):
#     print(f"{motif}: {result}")
