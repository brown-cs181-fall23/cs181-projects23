from hmmlearn import hmm
import numpy as np
import pickle


def train_hmm(motif, motif_pickle):
    """
    Train HMM on motif and negative control sequences and return trained HMM.
    Input: motif (str), motif_pickle (str)
    Output: hmm.CategoricalHMM
    """

    # Import Pickle File Here
    with open(motif_pickle, "rb") as f:
        motif_dict = pickle.load(f)

    # TODO: Set up HMM
    states = ["background", "motif"]

    # Do NOT change the random state. There will be other parameters to this function.
    model = hmm.CategoricalHMM(random_state=1)

    # TODO: Prepare Data for Training

    # TODO: Fit HMM

    return model


def predict(sequence) -> bool:
    """
    The base predict function. Predicts true if any predicted state is "motif".
    Feel free to change this!
    """
    return True if 1 in sequence else False


def test_hmm(model, data, ground_truths=None):
    """
    Given a set of sequences, returns whether they are likely to contain the motif or only background.
    Input: hmm.CategoricalHMM, sequences (list[string]), labels (list[int])
    Output: float
    """

    # TODO: Prepare Data for Decoding

    # TODO: Extract most likely sequence of states for each observation
    _, states_sequence = model.decode()

    # TODO: Currently the States Are a 1D array...reshape!

    # TODO: For each sequence, use the predict function to make a list of 1s and 0s.
    # 1: sequence containing a motif, 0: sequence missing a motif.
    test_results = []

    # TODO: Determine Model Performance (hint: np.sum())
    test_results = np.asarray(test_results)
    if ground_truths is not None:
        # Percent of obersvations that match ground truth labels
        pass
    else:
        # Check what fraction of sequences yielded positive predictions (label is 1)
        pass


def test_hmms(models, test_pickle, labeled):
    """
    Given multiple models for different motifs, tests a set of sequences and returns which motifs are likely to be present.
    Input: models (list[hmm.CategoricalHMM]), test_pickle (str), labeled (bool)
    Output: list[float]
    """
    with open(test_pickle, "rb") as f:
        motif_dict = pickle.load(f)

    test_results = []
    for motif in models:
        if labeled:
            observations, ground_truths = motif_dict[motif]
        else:
            # ground_truths will be None here
            observations, ground_truths = motif_dict["mystery"]
        test_results.append(test_hmm(models[motif], observations, ground_truths))
    return test_results
