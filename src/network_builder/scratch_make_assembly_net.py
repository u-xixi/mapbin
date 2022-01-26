import gfapy

# gfaf = "/Users/xixichen/Documents/projects/2019_haplotagging/02_assembly/sample1/assembly_graph_with_scaffolds.gfa"
# changed_gfaf = "/Users/xixichen/Documents/projects/2019_haplotagging/netbin/testing/gfa_nopaths.gfa"
# with open(gfaf, 'r') as f:
#     gcontent = iter(f.read().strip().split("\n"))
#
#
# new_out = []
# while True:
#
#     try:
#         curline = next(gcontent)
#         if not curline.startswith("P"):
#             new_out.append(curline)
#     except:
#         break
#
#
# with open(changed_gfaf, "w") as f:
#     f.write("\n".join(new_out) + "\n")

# changed_gfaf = "/Users/xixichen/Documents/projects/2019_haplotagging/netbin/testing/gfa_nopaths.gfa"
#
# with open(changed_gfaf, 'r') as f:
#     gcontent = f.read().strip()
# gfa_obj = gfapy.Gfa(gcontent)


from infomap import Infomap
for k, v in Infomap.__dict__.items():
    print(k)
