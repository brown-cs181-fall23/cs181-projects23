import hmmlearn as hmm
import numpy as np
from preprocess import preprocess_motif
from hmm_traintest import train_hmm, test_hmms
import pickle

## Train HMMs
motifs = ["ATATAT", "TATATA"]
bed_files = ["ATATAT.bed", "TATATA.bed"]
neg_file = "neg.bed"

## Read in Pickle Files
preprocess_motif(motifs, bed_files, neg_file)
with open("motif_dict.pkl", "rb") as f:
    motif_dict = pickle.load(f)

## Train HMM
trained_models = []
states_sequences = []
logprobs = []
for motif in motifs:
    model, logprob, states_sequence = train_hmm(motif, motif_dict)
    trained_models.append(model)
    states_sequences.append(states_sequence)
    logprobs.append(logprob)

## Test HMM
test_pickle = "test_dict.pkl"
with open(test_pickle, "rb") as f:
    test_dict = pickle.load(f)

results = test_hmms(trained_models, test_dict)

## Save Results

## Plot Results

## Save Plot
