import treeswift
import sys

t = treeswift.read_tree_newick(sys.argv[1])
rootcount = len(t.root.children)
for n in t.traverse_postorder():
	if n.is_root():
		continue
	if n.parent.is_root():
		n.edge_length=rootcount
		rootcount -= 1
	else:
		n.edge_length=-int(n.label)+int(n.parent.label)
print(t)
