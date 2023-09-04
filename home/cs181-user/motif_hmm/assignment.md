## Explore Informatics: Training Hidden Markov Models

Students will train a HMM to recognize sequences that contain a given motif. This will be done using the `hmmlearn` python package. 

There will be several sets of data that students will be provided that associate with a certain motif. Once student's have implemented the code to train the HMM, not write the actual training logic itself, they can create models that detect each motif given the provided labeled sequence and sequences that are associated with background regions. 

Then, the students will be provided with a random of sample of sequences containing a specific set of motifs. There will be a specific set of
motifs that together will suggest some function that students will have to look into. (This part can be workshopped).


Project Contents:
* There are three main files: `preprocess.py`, `hmm_traintest.py`, `explore_hmm.py`. These files are all appropriately commented.
* During TA Camp we can decide which functions and regions to abstract for students to figure out on their own. 
* The non-chat-GPTable portion of this projects should have students visualize the states and see whether the regions where motifs are detected actually reflect the length of the motif (i.e. is state 1 - motif present - for the length of the motif)
* Finally, students should reflect on the downfalls of this approach and errors and for extra credit can implement an alternate approach that achieves better motif recognition performance. 

Key Points:
In noisy data, it is possible that a motif is not actually detected because it is not an exact match. By training an HMM, we can detect whether a sequence may contain a motif even if the motif is not actually entirely there. 