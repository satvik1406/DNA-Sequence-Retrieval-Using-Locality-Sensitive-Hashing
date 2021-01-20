import numpy as np
from lsh import *


def lsh_query(signatures, num_bands):
    """
    Function for running the Locality Sensitive Hashing algorithm on the query as entered by user.
    (INPUT : Signatures as a list of lists , an integer.)
    (OUTPUT :a set of buckets.)
    """
    split_signatures = np.array([np.array_split(row, num_bands) for row in signatures])
    candidates = set()
    for j in range(num_bands):
        bucket = tuple(split_signatures[-1][j])
        for i in range(len(split_signatures)-1):
            if tuple(split_signatures[i][j]) == bucket:
                candidates.add((i,len(split_signatures)-1))
    return candidates

shingle_length = 6
num_hashes = 100
num_bands = 5
threshold = 0.9


"""
    Code to call various functions, determine time taken by each function and obtain results from given data
"""
while True:
    query = input("\n\nEnter a search query (Enter 'exit' to exit program):\t")
    if query == "exit":
        break

    start = time.time()
    t = time.time()
    print("Loading dataset...")
    data = load_dataset()
    print("Done.", len(data), "sequences loaded in", time.time()-t, "seconds.\n")

    t = time.time()
    print("Making shingles...")
    shingled_data = make_shingles(data, shingle_length)
    query_shingles = string_to_shingles(query, shingle_length)
    shingled_data.append(query_shingles)    # Append query to end of data
    print("Done.", "Took", time.time()-t, "seconds.\n")

    t = time.time()
    print("Minhashing...")
    signatures = minhash(shingled_data, num_hashes)
    print("Done.", "Signatures for", num_hashes, "hash functions generated in", time.time()-t, "seconds.\n")

    t = time.time()
    print("Generating candidates using LSH algorithm...")
    candidates = lsh_query(signatures, num_bands)
    print("Done.", len(candidates), "candidates generated in", time.time()-t, "seconds.\n")

    t = time.time()
    print("Checking similarity for candidates...")
    similar_strings = check_candidates(signatures, candidates, threshold)
    print("Done.", len(similar_strings), "similar sequences found in", time.time()-t, "seconds.\n")

    print("Total time taken: ", time.time()-start, "seconds.")
    if len(candidates)>0:
        precision = float(len(similar_strings))/len(candidates)
    else:
        precision = 0
    print("Precision: ", precision)
    print("\nShowing resultant similar documents:")
    for pair in similar_strings:
        # print(pair[0], ". ", data[pair[0]][0])
        print(pair[0])
