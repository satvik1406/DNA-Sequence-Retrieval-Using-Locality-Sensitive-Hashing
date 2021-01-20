import binascii
import numpy as np
import time
import random
import itertools
import matplotlib.pyplot as plt


def load_dataset():
    """
    fuction to load clethe dataset.
    The dataset has information regarding DNA secquence and the class of the DNA sequence.
    (OUTPUT : all the DNA sequences as a numpy array.)
    """
    data = []
    with open('./dataset/human_data.txt', 'r', encoding="utf-8") as f:
            out = f.readlines()
            for line in out:
                line = line.strip('\n').split()
                data.append([line[0], int(line[1])])
    return np.array(data)


def string_to_shingles(string, k):
    """
    Function to convert the document into k-shingles.
    (INPUT : string of charecters , shingle length as integer.)
    (OUTPUT : a list of unique hashed k-shingles of the input string.)
    """
    string = string.ljust(k)
    shingles = {string[i:i + k] for i in range(len(string) - k + 1)}
    shingles_hash = {binascii.crc32(bytes(shingle,'utf-8')) & 0xffffffff for shingle in shingles}
    return list(shingles_hash)


def make_shingles(data, k):
    """
    Function to convert the DNA sequence input into k-shingles.
    (INPUT : a tuple having DNA secquence and the class of DNA secquence , shingle length as integer.)
    (OUTPUT : a list of lists containing unique hashed k-shingles.)
    """
    return [string_to_shingles(row[0], k) for row in data]

def generate_hash_functions(num_hashes):
    """
    Function to genrate the hash functions.
    (INPUT : total number of hash functions as an integer.)
    (OUTPUT : 2 lists of random integers and a large prime number which constitute to our hash functions.)
    """
    a = set()
    b = set()
    while len(a)<num_hashes:
        a.add(random.randint(0, 2**32 - 1))
    while len(b)<num_hashes:
        b.add(random.randint(0, 2**32 - 1))
    c = 4294967311
    return list(a), list(b), c

def minhash(shingled_data, num_hashes):
    """
    Function for minhashing the shingled data.
    (INPUT : shingled data , the total number of hash functions.)
    (OUTPUT : Signatures as a list of lists.)
    """
    a, b, c = generate_hash_functions(num_hashes)
    signatures = []
    for shingle_set in shingled_data:
        signature = []
        for i in range(num_hashes):
            shingle_set = np.array(shingle_set)
            signature.append(np.min((shingle_set*a[i] + b[i])%c))
        signatures.append(signature)
    return signatures

def lsh(signatures, num_bands):
    """
    Function for running the Locality Sensitive Hashing algorithm.
    (INPUT : Signatures as a list of lists , an integer.)
    (OUTPUT :a set of buckets.)
    """

    split_signatures = np.array([np.array_split(row, num_bands) for row in signatures])
    candidates = set()
    for j in range(num_bands):
        buckets = {}
        for i in range(len(split_signatures)):
            if tuple(split_signatures[i][j]) not in buckets:
                buckets[tuple(split_signatures[i][j])] = []
            buckets[tuple(split_signatures[i][j])].append(i)
        for row in buckets.values():
            candidates |= set(list(itertools.combinations(row,2)))
        buckets.clear()
    return candidates

def hamming_similarity(signature_1, signature_2):
    """
    Function for finding the hamming similarity betweeen two signatures.
    (INPUT : two signature lists.)
    (OUTPUT : hamming similarity of the inputs.)
    """
    count = 0
    for a,b in zip(signature_1, signature_2):
        if a==b:
            count+=1
    return float(count)/float(len(signature_1))

def check_candidates(signatures, candidates, threshold):
    """
    Function to collect pairs of similar sequences.
    (INPUT : signatures , set of buckets , the threshold value for similarity comparision.)
    (OUTPUT : all similar pairs as list.)
    """
    similar_pairs = []
    for pair in candidates:
        if pair[0] == pair[1]:
            flag += 1
        if hamming_similarity(signatures[pair[0]], signatures[pair[1]]) > threshold:
            similar_pairs.append(pair)
    return similar_pairs

def make_plots():
    """
    Function to plot graphs to determine the optimal values of number of bands , shingle length(k) etc.
    """
    # To print precision graphs and graph of number of candidate vs similar pairs
    x = []
    y1 = []
    y2 = []
    for threshold in range(50,100,5):
        shingle_length = 10
        num_hashes = 100
        num_bands = 5
        threshold = float(threshold)/100

        data = load_dataset()
        shingled_data = make_shingles(data, shingle_length)
        signatures = minhash(shingled_data, num_hashes)
        candidates = lsh(signatures, num_bands)
        similar_pairs = check_candidates(signatures, candidates, threshold)
        x.append(threshold)
        y1.append(len(candidates))
        y2.append(len(similar_pairs))

    y = np.divide(y2,y1)

    plt.figure()
    plt.plot(x, y)
    plt.xlabel('Threshold')
    plt.ylabel('Precision')
    plt.legend()
    plt.savefig("threshold_vs_precision")
    plt.close()

    plt.figure()
    plt.plot(x, y1, label = "Candidate pairs")
    plt.plot(x, y2, label = "Similar pairs")
    plt.xlabel('Threshold')
    plt.ylabel('Number of pairs')
    plt.savefig("threshold_vs_pairs")
    plt.close()

    # To print graph of Similarity vs Probability of being detected as candidate
    plt.figure()
    fig, ax = plt.subplots(2,3)
    
    ax[0,0].set(xlabel='Similarity', ylabel='Probability')
    ax[0,0].set_title('b = [5,10...50], r = 5')
    r=5
    for b in range(5,55,5):
        x = []
        y = []
        for s in np.arange(0.01,1,0.01):
            prob = 1-((1-s**r)**b)
            x.append(s)
            y.append(prob)
        ax[0,0].plot(x, y, c = 'purple', lw = 0.5)

    ax[0,1].set(xlabel='Similarity', ylabel='Probability')
    ax[0,1].set_title('b = 20, r = [5,10...50]')
    b=50
    for r in range(5,55,5):
        x = []
        y = []
        for s in np.arange(0.01,1,0.01):
            prob = 1-((1-s**r)**b)
            x.append(s)
            y.append(prob)
        ax[0,1].plot(x, y, c = 'purple', lw = 0.5)

    ax[0,2].set(xlabel='Similarity', ylabel='Probability')
    ax[0,2].set_title('h = 100, b = [5,10...50], r = h/b')
    h=100
    for b in range(5,50,5):
        r = h/b
        x = []
        y = []
        for s in np.arange(0.01,1,0.01):
            prob = 1-((1-s**r)**b)
            x.append(s)
            y.append(prob)
        ax[0,2].plot(x, y, c = 'purple', lw = 0.5)

    ax[1,0].set(xlabel='Similarity', ylabel='Probability')
    ax[1,0].set_title('h = 200, b = [5,10...100], r = h/b')
    h=200
    for b in range(5,100,5):
        r = h/b
        x = []
        y = []
        for s in np.arange(0.01,1,0.01):
            prob = 1-((1-s**r)**b)
            x.append(s)
            y.append(prob)
        ax[1,0].plot(x, y, c = 'purple', lw = 0.5)

    ax[1,1].set(xlabel='Similarity', ylabel='Probability')
    ax[1,1].set_title('h = 400, b = [10,20...200], r = h/b')
    h=400
    for b in range(10,200,10):
        r = h/b
        x = []
        y = []
        for s in np.arange(0.01,1,0.01):
            prob = 1-((1-s**r)**b)
            x.append(s)
            y.append(prob)
        ax[1,1].plot(x, y, c = 'purple', lw = 0.5)

    ax[1,2].set(xlabel='Similarity', ylabel='Probability')
    ax[1,2].set_title('h = 500, b = [10,20...200], r = h/b')
    h=500
    for b in range(10,200,10):
        r = h/b
        x = []
        y = []
        for s in np.arange(0.01,1,0.01):
            prob = 1-((1-s**r)**b)
            x.append(s)
            y.append(prob)
        ax[1,2].plot(x, y, c = 'purple', lw = 0.5)
        

    plt.tight_layout()
    plt.show()
    plt.close()

# make_plots()
