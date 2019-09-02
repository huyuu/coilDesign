# Gauss Integral Nodes and Weights
import Pkg
Pkg.add("FastGaussQuadrature")
using FastGaussQuadrature
nodes, weights = gausslegendre(100)
println(nodes)
println(weights)

fileNodes = open("gaussNodesForPrecise.csv", "w")
fileWeights = open("gaussWeightsForPrecise.csv", "w")

for (node, weight) in zip(nodes, weights)
    write(fileNodes, "$node\n")
    write(fileWeights, "$weight\n")
end

close(fileNodes)
close(fileWeights)
