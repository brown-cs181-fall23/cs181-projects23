import hmmlearn as hmm
import numpy as np
import pickle


def train_hmm(motif, motif_pickle):
    """
    Train HMM on motif and negative control sequences and return trained HMM.
    Input: motif (str), motif_pickle (str)
    Output: hmm.MultinomialHMM
    """
    ## Import Pickle File Here
    with open("motif_dict.pkl", "rb") as f:
        motif_dict = pickle.load(f)
    observations = motif_dict[motif]

    ## Train HMM
    ## Initialize HMM
    states = ["background", "motif"]
    id2topic = dict(enumerate(states))
    start_probs = np.array([0.9, 0.1])

    vocabulary = ["A", "T", "C", "G"]
    # Convert nucleotides to IDs
    nucleotide2id = dict(zip(vocabulary, range(len(vocabulary))))

    # Convert sequences to numerical representation
    X = []
    for sequence in observations:
        row = [nucleotide2id[nucleotide] for nucleotide in sequence]
        X.append(row)

    data = np.array(X, dtype=int)

    # Set up model
    n_states = len(states)
    n_nucleotides = len(vocabulary)
    n_trials = len(observations[0]) * len(observations)  # Length of each sequence

    model = hmm.MultinomialHMM(
        n_components=n_states, n_trials=n_trials, n_iter=50, init_params=""
    )

    model.n_features = n_nucleotides
    model.startprob_ = np.array([0.5, 0.5])  # Adjust as needed
    model.transmat_ = np.array([[0.8, 0.2], [0.2, 0.8]])  # Adjust as needed

    # Initialize emission probabilities
    model.emissionprob_ = np.array(
        [[0.25, 0.25, 0.25, 0.25], [0.25, 0.25, 0.25, 0.25]]  # motif
    )  # background

    model.fit(data, [n_trials] * len(observations))
    logprob, states_sequence = model.decode(data)
    return model, logprob, states_sequence


def test_motif(model, data):
    """
    Given a set of sequences, return whether they are likely to contain the motif or only background.
    Input: hmm.MultinomialHMM, sequences (list)
    Output: list : Boolean
    """
    ## Test HMM
    logprob, states_sequence = model.decode(data)

    ## For each sequence in sequences, return whether it is likely to contain the motif or only background
    test_results = []
    for sequence_states in states_sequence:
        if 1 in sequence_states:
            test_results.append(True)
        else:
            test_results.append(False)
    return test_results


def test_hmms(models, sequences):
    """
    Given multiple models for different motifs, test a set of sequences and return which motifs are likely to be present.
    """
    ## Convert sequences to numerical representation
    vocabulary = ["A", "T", "C", "G"]
    # Convert nucleotides to IDs
    nucleotide2id = dict(zip(vocabulary, range(len(vocabulary))))

    X = []
    for sequence in sequences:
        row = [nucleotide2id[nucleotide] for nucleotide in sequence]
        X.append(row)

    data = np.array(X, dtype=int)

    ## For each model, test sequences
    test_results = []
    for model in models:
        test_results.append(test_motif(model, sequences, data))
    return test_results
