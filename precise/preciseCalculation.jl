# This is the precise prediction for coil design.
using Distributed
addprocs(4)
@everywhere using Statistics

# Constants
# Physical Constants
@everywhere const μ = 4pi*1e-7
# Current
@everywhere const I = 3  # 1[A]
@everywhere const N = 400
# Ingredient
@everywhere const standardPhiOfConductor = 1.02e-3  # 1.02mm using AWG 20
@everywhere const thicknessOfGFRPWall = 2e-2  # 2cm
# Measurement Area
@everywhere const X0 = 1e-2  # 1cm
@everywhere const Y0 = 1e-2  # 1cm
@everywhere const Z0 = 1e-2  # 1cm
# Coil Shape
@everywhere const h = 5e-2 # 5cm
@everywhere const l = 2h # 5cm
@everywhere const R = h # 5cm
@everywhere const d = 0.5h # 5cm


# Variables

# Measurement points
@everywhere const sampleIntervals = 20
@everywhere const samplePoints = sampleIntervals+1
# Gauss Integral Nodes and Weights
@everywhere const nodes = let
    nodes = []
    open("gaussNodesForPrecise.csv", "r") do file
        lines = readlines(file)
        for line in lines
            newNode = parse(Float64, line)
            push!(nodes, newNode)
        end
    end
    nodes
end
@everywhere const weights = let
    weights = []
    open("gaussNodesForPrecise.csv", "r") do file
        lines = readlines(file)
        for line in lines
            newWeight = parse(Float64, line)
            push!(weights, newWeight)
        end
    end
    weights
end


# Children Variables
@everywhere const standardRadiusOfConductor = standardPhiOfConductor/2
@everywhere const conductorsPerLayer = div(thicknessOfGFRPWall, standardPhiOfConductor)
# sample points
@everywhere const xs = LinRange(-X0, X0, samplePoints)
@everywhere const ys = LinRange(-Y0, Y0, samplePoints)
@everywhere const zs = LinRange(-Z0, Z0, samplePoints)
# File Operation
const dirName = "I=$(round(I, sigdigits=2))_N=$(round(Int, N))_h=$(round(h*100, sigdigits=2))cm_l=$(round(l*100, sigdigits=2))cm_X0=$(round(X0*100, sigdigits=2))cm_Y0=$(round(Y0*100, sigdigits=2))cm_Z0=$(round(Z0*100, sigdigits=2))cm_conductorPhi=$(round(standardPhiOfConductor*1e3, sigdigits=3))mm"


# Models

@everywhere const BVector = Array{Float64, 1}

# Shape parameters of source elements.
@everywhere struct ParsOfSource
    c1_y::Float64
    c1_z::Float64
    c2_y::Float64
    c2_z::Float64
    c3_R::Float64
    c3_z::Float64
    c4_R::Float64
    c4_z::Float64

    c5_y::Float64
    c5_z::Float64
    c6_y::Float64
    c6_z::Float64
    c7_R::Float64
    c7_z::Float64
    c8_R::Float64
    c8_z::Float64
end
# Generate from turn number n.
@everywhere ParsOfSource(n::Int, h::Float64) = let
    index::PositionIndex = PositionIndex(n)

	global d
    c1_y =  -h - standardPhiOfConductor*(index.layer-1) - standardRadiusOfConductor
    c1_z = -d + index.centerDistance * standardPhiOfConductor
    c2_y = h + standardPhiOfConductor*(index.layer-1) + standardRadiusOfConductor
    c2_z = c1_z

    R = h
    c3_R = R + standardPhiOfConductor*(index.layer-1) + standardRadiusOfConductor
    c3_z = c1_z
    c4_R = c3_R
    c4_z = c1_z

    c5_y =  -h - standardPhiOfConductor*(index.layer-1) - standardRadiusOfConductor
    c5_z = d + index.centerDistance * standardPhiOfConductor
    c6_y = h + standardPhiOfConductor*(index.layer-1) + standardRadiusOfConductor
    c6_z = c5_z

    c7_R = R + standardPhiOfConductor*(index.layer-1) + standardRadiusOfConductor
    c7_z = c5_z
    c8_R = c7_R
    c8_z = c5_z

    return ParsOfSource(c1_y, c1_z, c2_y, c2_z, c3_R, c3_z, c4_R, c4_z, c5_y, c5_z, c6_y, c6_z, c7_R, c7_z, c8_R, c8_z)
end


@everywhere struct PositionIndex
    layer::Int
    item::Int
    centerDistance::Float64
end
@everywhere PositionIndex(n::Int) = let
    layer = cld(n, conductorsPerLayer)
    item = n % conductorsPerLayer
    center = (conductorsPerLayer+1)/2
    centerDistance::Float64 = item - center
    PositionIndex(layer, item, centerDistance)
end


struct ResultsOfMeanVarRate
    meanBVector::BVector
    varRateVector::BVector
end
ResultsOfMeanVarRate(meanBVector::BVector, minBVector::BVector, maxBVector::BVector) = let
    varRateVector = (maxBVector - minBVector) / meanBVector
    ResultsOfMeanVarRate(meanBVector, varRateVector)
end


@everywhere function numericalIntegrateOf(f::Function)::BVector
    local result::BVector = [0, 0, 0]
    for (node, weight) in zip(nodes, weights)
        result += weight * f(node)
    end
    return result
end


@everywhere function c1(x::Float64, y::Float64, z::Float64; u::Float64, _y::Float64, _z::Float64, l::Float64)::BVector
    xVector::Float64 = 0
    yVector::Float64 = -(z-_z)
    zVector::Float64 = y - _y
    denominator::Float64 = ( (x-l*u)^2 + (y-_y)^2 + (z-_z)^2 ) ^ (1.5)
    return [xVector, l*yVector/denominator, l*zVector/denominator]
end


@everywhere function c2(x::Float64, y::Float64, z::Float64; u::Float64, _y::Float64, _z::Float64, l::Float64)::BVector
    xVector::Float64 = 0
    yVector::Float64 = z-_z
    zVector::Float64 = -(y - _y)
    denominator::Float64 = ( (x-l*u)^2 + (y-_y)^2 + (z-_z)^2 ) ^ (1.5)
    return [xVector, l*yVector/denominator, l*zVector/denominator]
end


@everywhere function c3(x::Float64, y::Float64, z::Float64; R::Float64, u::Float64, _z::Float64, l::Float64)::BVector
    phi::Float64 = 0.5pi*u
    xVector::Float64 = R*cos(phi)*(z-_z)
    yVector::Float64 = R*sin(phi)*(z-_z)
    zVector::Float64 = -R*( (y-R*sin(phi))sin(phi) + (x-l-R*cos(phi))cos(phi) )
    denominator::Float64 = ( (x-l-R*cos(phi))^2 + (y-R*sin(phi))^2 + (z-_z)^2 ) ^ (1.5)
    return 0.5pi * [xVector/denominator, yVector/denominator, zVector/denominator]
end


@everywhere function c4(x::Float64, y::Float64, z::Float64; R::Float64, u::Float64, _z::Float64, l::Float64)::BVector
    phi::Float64 = 0.5pi*(u+2)
    xVector::Float64 = R*cos(phi)*(z-_z)
    yVector::Float64 = R*sin(phi)*(z-_z)
    zVector::Float64 = -R*( (y-R*sin(phi))sin(phi) + (x+l-R*cos(phi))cos(phi) )
    denominator::Float64 = ( (x+l-R*cos(phi))^2 + (y-R*sin(phi))^2 + (z-_z)^2 ) ^ (1.5)
    return 0.5pi * [xVector/denominator, yVector/denominator, zVector/denominator]
end


@everywhere function BFromSingleTurn(x, y, z; n::Int, h::Float64, l::Float64)::BVector
    pars::ParsOfSource = ParsOfSource(n, h)
    local result::BVector = [0, 0, 0]
    # B from lower
    result .+= numericalIntegrateOf((u) -> c1(x, y, z; u=u, _y=pars.c1_y, _z=pars.c1_z, l=l))
    result .+= numericalIntegrateOf((u) -> c2(x, y, z; u=u, _y=pars.c1_y, _z=pars.c1_z, l=l))
    result .+= numericalIntegrateOf((u) -> c3(x, y, z; R=pars.c3_R, u=u, _z=pars.c3_z, l=l))
    result .+= numericalIntegrateOf((u) -> c4(x, y, z; R=pars.c4_R, u=u, _z=pars.c4_z, l=l))
    # B from upper
    result .+= numericalIntegrateOf((u) -> c1(x, y, z; u=u, _y=pars.c5_y, _z=pars.c5_z, l=l))
    result .+= numericalIntegrateOf((u) -> c2(x, y, z; u=u, _y=pars.c6_y, _z=pars.c6_z, l=l))
    result .+= numericalIntegrateOf((u) -> c3(x, y, z; R=pars.c7_R, u=u, _z=pars.c7_z, l=l))
    result .+= numericalIntegrateOf((u) -> c4(x, y, z; R=pars.c8_R, u=u, _z=pars.c8_z, l=l))
    return result
end


@everywhere function calculateBAtPoint(x::Float64, y::Float64, z::Float64; h::Float64, l::Float64)::BVector
    local results::BVector = zeros(3)

    map(1:N) do n
        results .+= BFromSingleTurn(x, y, z; n=n, h=h, l=l)
    end

    return μ/(4pi)*I * results
end


function calculateResultWhen(; h::Float64, l::Float64)::ResultsOfMeanVarRate
	# init results
    meanBVector::BVector = zeros(3)
    minBVector::BVector = ones(3)
    maxBVector::BVector = zeros(3)
	# for calculation
    totalSamplePoints = length(xs)*length(ys)*length(zs)
    futureBVectorsInCube::Array{Future} = []

    for (xIndex, xValue) in enumerate(xs), (yIndex, yValue) in enumerate(ys), (zIndex, zValue) in enumerate(zs)
        future = @spawn calculateBAtPoint(xValue, yValue, zValue; h=h, l=l)
        push!(futureBVectorsInCube, future)
    end

    map(futureBVectorsInCube) do futureBVector
        bVector::BVector = fetch(futureBVector)
        bVector = abs.(bVector)
        meanBVector += bVector
        map(1:3) do i
            minBVector[i] = bVector[i] < minBVector[i] ? bVector[i] : minBVector[i]
            maxBVector[i] = bVector[i] > maxBVector[i] ? bVector[i] : maxBVector[i]
        end
    end
    meanBVector /= totalSamplePoints

	println("maxBVector = $(maxBVector)")
	println("minBVector = $(minBVector)")
    return ResultsOfMeanVarRate(meanBVector, minBVector, maxBVector)
end


function openWithNewDirIfNeeded(;dirName::String, fileName::String, modes::String)
    isSuccessful = true
    file = try
        open("$dirName/$fileName", modes)
    catch
        mkdir(dirName)
        touch("$dirName/$fileName")
        isSuccessful = false
    end

    if isSuccessful == false
        file = open("$dirName/$fileName", modes)
    end
    return file
end


function myOpen(;fileName::String, modes::String="w", dirName::Union{String, Nothing}=nothing, withNewDirIfNeeded::Bool=true, csvHeader::Union{String, Nothing}=nothing)
    if withNewDirIfNeeded
        if isa(dirName, String)
            file = openWithNewDirIfNeeded(;dirName=dirName, fileName=fileName, modes=modes)
        else  # dirName == nothing
            error("dirName shouldn't be nothing.")
        end
    end

    if isa(csvHeader, String)
        write(file, "$(csvHeader)\n")
    end

    return file
end


# Main

let
    result = @time calculateResultWhen(; h=h, l=l)
    resultFile = myOpen(;fileName="resultsWhenSamples=$(sampleIntervals).csv", modes="w", dirName=dirName, csvHeader="varRateX[%],varRateY[%],varRateZ[%],meanBx[mT],meanBy[mT],meanBz[mT]")
    write(resultFile, "$(result.varRateVector[1]*100),")
    write(resultFile, "$(result.varRateVector[2]*100),")
    write(resultFile, "$(result.varRateVector[3]*100),")
    write(resultFile, "$(result.meanBVector[1]*1000),")
    write(resultFile, "$(result.meanBVector[2]*1000),")
    write(resultFile, "$(result.meanBVector[3]*1000)\n")
    close(resultFile)
end
