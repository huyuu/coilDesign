# Gauss Integral Nodes and Weights
import Pkg
Pkg.add("FastGaussQuadrature")
using FastGaussQuadrature
nodes, weights = gausslegendre(10)
println(nodes)
println(weights)

fileNodes = open("gaussNodes.csv", "w")
fileWeights = open("gaussWeights.csv", "w")

for (node, weight) in zip(nodes, weights)
    write(fileNodes, "$node\n")
    write(fileWeights, "$weight\n")
end

close(fileNodes)
close(fileWeights)
