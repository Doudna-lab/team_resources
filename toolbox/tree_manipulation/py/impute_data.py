from Bio import Phylo
import gzip
import pylab


def main():
	# with open(", 'rb') as tree:
	tree = "/Users/bellieny/projects/team_resources/toolbox/tree_manipulation/metadata/ar53_r214.tree"
	nwk = Phylo.read(tree, "newick")
	nwk.ladderize()
	Phylo.draw_graphviz(nwk)
	pylab.show()

	clades = []
	for i in nwk:
		clades.append(i)
		print(i)

if __name__ == "__main__":
	main()
