input_dna = "AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC"

a = 0
c = 0
t = 0
g = 0

for i in input_dna:
    if i == "A":
        a+=1
    elif i == "C":
        c +=1
    elif i == "T":
        t += 1
    elif i == "G":
        g += 1

print(a, c, t, g)