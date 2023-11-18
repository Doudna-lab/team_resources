from Bio import Phylo
import gzip
import pylab


def get_parent(tree, child_clade):
    node_path = tree.get_path(child_clade)
    return node_path[-2]


def main():
	# with open(", 'rb') as tree:
	tree = "/Users/bellieny/projects/team_resources/toolbox/tree_manipulation/metadata/ar53_r214.tree"
	nwk = Phylo.read(tree, "newick")
	nwk.ladderize()
	Phylo.draw_graphviz(nwk)
	pylab.show()

	# Select a clade
	myclade = nwk.find_clades(1)
	# Test the function
	parent = get_parent(nwk, myclade)
	assert myclade in parent

	clades = []
	for i in nwk:
		clades.append(i)
		print(i)

if __name__ == "__main__":
	main()
