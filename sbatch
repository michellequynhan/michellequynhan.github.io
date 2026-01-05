from collections import defaultdict
import glob
import subprocess
import re
import sys
import os

treeDict = {}
bicdict = {}
failedlogdict = {}

# creating a dictionary of [geneID][speciesname]:[genome]
def dictMaker(fileList):
    seqdict = defaultdict(dict)
    for x in fileList:
        file1 = open(x, 'r')
        speciesname = x.split('.')[0]
        for line in file1:
            z = len(line)
            if line[0:z-1] == '>' + geneID:
                seqdict[geneID][speciesname] = file1.readline().replace('\n', '')
    return seqdict

# filtering out all dash (-) or X genomes because they have no information and break my code
def dash_finder(dictA):
    values = dictA[geneID].values()
    comp_list = []
    pattern = r'[^\-XN]'
    for x in values:
        if re.search(pattern, x):
            comp_list.append(False)
        else:
            comp_list.append(True)
        all_list = any(comp_list)
    return all_list

# filtering out genes that have identical genomes across all species
def equalizer(dictA):
    values = dictA[geneID].values()
    lis = list(values)
    res = False
    if (lis.count(lis[0]) + 1 == len(lis)) or (lis.count(lis[0]) == len(lis)):
        res = True
        return(str(res))
    else:
        print('False')

dict_list = ['genus', 'names', 'censored', 'for', 'data', 'privacy']

# creating the multiple sequence alignment from the first dictionary
def write_msa(dictA):
    f = open(geneID + '.msa.txt', 'w+')
    for speciesname in dict_list:
        f.write('>' + speciesname + '\n')
        f.write(dictA[geneID][speciesname] + '\n')

# running iqtree
def run_iqtree(msafile):
    iqtree_cmd = 'iqtree -s ' + msafile + ' -B 1000'
    subprocess.call(iqtree_cmd, shell = True)

# running iqtree with a tree to be compared to
def two_iqtree(msafile):
    iqtree_cmd2 = 'iqtree -s ' + msafile + ' -m ' + bicdict[geneID] + ' -t ref_tree --prefix ' + geneID + '.top -redo'
    subprocess.call(iqtree_cmd2, shell = True)

# pulling the tree file for each gene (the text version of a phylogenetic tree)
def opentreefile(geneID):
    f = open(treefile, 'r')
    tree = f.read()
    return(tree)

# pulling the tree file (aligned to a specific example tree for a specific topology)
def opentoptreefile(geneID):
    f = open(toptreefile, 'r')
    tree = f.read()
    return(tree)

# pulling important info from iqtree output file
def openlogbic(geneID):
    f = open(logfile, "r")
    BIC = re.search("(?<=Bayesian Information Criterion:).*", f.read())
    print(BIC)
    BIC1 = str(BIC)
    BIC2 = BIC1.split("'")[1]
    BIC3 = BIC2.strip()
    return(BIC3)

# pulling failed genes
def failed_openlog(geneID):
    f = open(logfile, "r")
    species = re.findall(".*\%    failed", f.read())
    for i, s in enumerate(species):
        if s[2] == " ":
            species[i] = s.split(" ")[5]
        else:
            species[i] = s.split(" ")[4]
    return(species)

# creating gene list for forloop to run through
def geneListmk(file):
    geneListNew = []

    a = open(file, 'r')
    for line in a:
        if line[0] == '>':
            name = line.strip()
            name1 = name[1:]
            geneListNew.append(name1)
    return geneListNew

# changing to my own temporary directory in cluster so as to not complete these resource intensive jobs on the main cluster that everyone in the lab shares
tmp_cmd = 'mkdir /tmp/fish'
subprocess.call(tmp_cmd, shell = True)

# since there were so many genes, I split all files into 20 smaller files to be ran at the same time in parallel, rather than one giant job. see below for splitter script
string_arg = str(sys.argv[1]).strip()

# copying everything i need to the temporary directory
cp_cmd = 'cp *.fa' + string_arg + ' /tmp/fish'
subprocess.call(cp_cmd, shell = True)

noto_cp = 'cp ref_tree /tmp/fish'
subprocess.call(noto_cp, shell = True)

os.chdir('/tmp/fish')

failed_gene = []
w = open('treefile' + string_arg + '.data', 'w+')
b = open('bic_failed' + string_arg + '.data', 'w+')
g = open('failed_species' + string_arg + '.data', 'w+')
t = open('top_treefile' + string_arg + '.data', 'w+')

genus1new = 'genus1.fa' + string_arg
geneListNew = geneListmk(aceratusnew)

for geneID in geneListNew:
    #compiling all files after they have been split
    fileList = glob.glob('*.fa'  + string_arg)
    dictA = dictMaker(fileList)
    #filtering for genes that have dashes, Ns, Xs, or anything else not useful
    any_list = dash_finder(dictA)
    if any_list == True:
        print('Skipped '+ geneID)
        failed_gene.append(geneID)
    else:
        #good genes move to the next filter, which checks if all the genomes across species are the same for that one gene (there can't be convergence if every base is the same!)
        if equalizer(dictA) == 'True':
            print('Skipped ' + geneID)
            failed_gene.append(geneID)
        else:
            #genes that passed both tests run through iqtree twice
            write_msa(dictA)
            msafile = geneID + '.msa.txt'

            run_iqtree(msafile)
            #pulling the tree and other important information from the iqtree runs
            treefile = geneID + '.msa.txt.treefile'
            gene_tree = opentreefile(treefile)
            treeDict[geneID] = gene_tree
            w.write('>' + geneID + '\n')
            w.write(gene_tree)

            logfile = geneID + '.msa.txt.log'
            return_bic = openlogbic(logfile)
            bicdict[geneID] = return_bic
            return_failed = failed_openlog(logfile)
            failedlogdict[geneID] = return_failed
            failed_list = ','.join(return_failed)
            geneIDs = geneID.ljust(20)
            return_bics = return_bic.ljust(10)
            b.write(geneIDs + '\t' + return_bics + '\t' + 'Failed species: ' + failed_list + '\n')
            two_iqtree(msafile)
            toptreefile = geneID + '.top.treefile'
            top_gene_tree = opentoptreefile(toptreefile)
            t.write('>' + geneID + '\n')
            t.write(top_gene_tree)
            # removing all unneeded files after processing to reduce strain on cluster
            remove_cmd = 'rm ' + geneID + '*'
            subprocess.call(remove_cmd, shell = True)


failed_str = str(failed_gene)
g.write('Failed genes: ' + failed_str)
# this part/file didn't work but it didn't matter for my results so I didn't bother trouble shooting

g.close()
w.close()
b.close()
t.close()

# copying results to main directory
cp_cmd = 'mv *' + string_arg + '.data /project/daane/michelle/concatenated'
subprocess.call(cp_cmd, shell = True)
