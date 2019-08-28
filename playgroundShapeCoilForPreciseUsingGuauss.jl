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
# Coil Shape
@everywhere const h = 0.05  # 5cm
@everywhere const R = h
@everywhere const d = 0.5h
@everywhere const l = 2h
# Measurement Area
@everywhere const X0 = 0.1h
@everywhere const Y0 = 0.1h
@everywhere const Z0 = 0.1h
# Gauss Integral Nodes and Weights
@everywhere import Pkg
@everywhere Pkg.add("FastGaussQuadrature")
@everywhere using FastGaussQuadrature
@everywhere const nodes, weights = gausslaguerre(10)


# Variables

# Intervals into which a current source loop is cut, ex: c1. Should not be too small otherwise divergence condition could be raised."
@everywhere const sourceIntervals = 100
# Measurement points
@everywhere const sampleIntervals = 100
@everywhere const samplePoints = sampleIntervals+1


# Children Variables

@everywhere const xs = LinRange(-X0, X0, samplePoints)
@everywhere const ys = LinRange(-Y0, Y0, samplePoints)
@everywhere const zs = LinRange(-Z0, Z0, samplePoints)
@everywhere const standardRadiusOfConductor = standardPhiOfConductor/2
@everywhere const conductorsPerLayer = div(thicknessOfGFRPWall, standardPhiOfConductor)
# File Operation
const dirName = "preciseGauss_I=$(round(I, sigdigits=2))_N=$(round(Int, N))_h=$(round(h*100, sigdigits=2))cm_X0=$(round(X0/h, sigdigits=2))h_Y0=$(round(Y0/h, sigdigits=2))h_Z0=$(round(Z0/h, sigdigits=2))h_thickness=$(round(thicknessOfGFRPWall*100))cm_conductorPhi=$(round(standardPhiOfConductor*1e3, sigdigits=3))mm"



# Models

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
@everywhere ParsOfSource(n::Int) = let
    index::PositionIndex = PositionIndex(n)

    c1_y =  -h - standardPhiOfConductor*(index.layer-1) - standardRadiusOfConductor
    c1_z = -d + index.centerDistance * standardPhiOfConductor
    c2_y = h + standardPhiOfConductor*(index.layer-1) + standardRadiusOfConductor
    c2_z = c1_z

    global R
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


@everywhere function numericalIntegrateOf(f::Function)::Float64
    # sumOfIntervals = let
    #     intervals = [ lowerLimit + k*(upperLimit-lowerLimit)/n for k in 1:n-1 ]
    #     sum( f(k) for k in intervals )
    # end
    # return (upperLimit-lowerLimit)/n * ( (f(upperLimit)+f(lowerLimit))/2 + sumOfIntervals )
    return sum( f.(nodes) .* f.(weights) )
end


@everywhere function c1(x::Float64, y::Float64, z::Float64; u::Float64, _y::Float64, _z::Float64)::Float64
    zVector = y - _y
    denominator = ( (x-l*u)^2 + (y-_y)^2 + (z-_z)^2 ) ^ (1.5)
    return l*zVector/denominator
end


@everywhere function c2(x::Float64, y::Float64, z::Float64; u::Float64, _y::Float64, _z::Float64)::Float64
    zVector = -(y - _y)
    denominator = ( (x-l*u)^2 + (y-_y)^2 + (z-_z)^2 ) ^ (1.5)
    return l*zVector/denominator
end


@everywhere function c3(x::Float64, y::Float64, z::Float64; R::Float64, u::Float64, _z::Float64)::Float64
    zVector = -R*( (y-R*sin(0.5pi*u))sin(0.5pi*u) + (x-l-R*cos(0.5pi*u))cos(0.5pi*u) )
    denominator = ( (x-l-R*cos(0.5pi*u))^2 + (y-R*sin(0.5pi*u))^2 + (z-_z)^2 ) ^ (1.5)
    return 0.5pi*zVector/denominator
end


@everywhere function c4(x::Float64, y::Float64, z::Float64; R::Float64, u::Float64, _z::Float64)::Float64
    zVector = -R*( (y-R*sin(0.5pi*(u+2)))sin(0.5pi*(u+2)) + (x+l-R*cos(0.5pi*(u+2)))cos(0.5pi*(u+2)) )
    denominator = ( (x+l-R*cos(0.5pi*(u+2)))^2 + (y-R*sin(0.5pi*(u+2)))^2 + (z-_z)^2 ) ^ (1.5)
    return 0.5pi*zVector/denominator
end


@everywhere function BFromSingleTurn(x, y, z; n::Int)::Float64
    pars::ParsOfSource = ParsOfSource(n)
    local result::Float64 = 0
    # B from lower
    result += numericalIntegrateOf((u) -> c1(x, y, z; u=u, _y=pars.c1_y, _z=pars.c1_z))
    result += numericalIntegrateOf((u) -> c2(x, y, z; u=u, _y=pars.c1_y, _z=pars.c1_z))
    result += numericalIntegrateOf((u) -> c3(x, y, z; R=pars.c3_R, u=u, _z=pars.c3_z))
    result += numericalIntegrateOf((u) -> c4(x, y, z; R=pars.c4_R, u=u, _z=pars.c4_z))
    # B from upper
    result += numericalIntegrateOf((u) -> c1(x, y, z; u=u, _y=pars.c5_y, _z=pars.c5_z))
    result += numericalIntegrateOf((u) -> c2(x, y, z; u=u, _y=pars.c6_y, _z=pars.c6_z))
    result += numericalIntegrateOf((u) -> c3(x, y, z; R=pars.c7_R, u=u, _z=pars.c7_z))
    result += numericalIntegrateOf((u) -> c4(x, y, z; R=pars.c8_R, u=u, _z=pars.c8_z))
    return result
end


function calculateBAt(x::Float64, y::Float64, z::Float64)::Float64
    local results::Float64 = 0
    futureResultsOfSingleTurn::Array{Future} = []

    map(1:N) do n
        future = @spawn BFromSingleTurn(x, y, z; n=n)
        push!(futureResultsOfSingleTurn, future)
    end

    map(futureResultsOfSingleTurn) do futureResult
        results += fetch(futureResult)
    end

    return μ/(4pi)*I * results
end


function openWithNewDirIfNeeded(;dirName::String, fileName::String, modes::String)
    isSuccessful = true
    file = try
        open("$dirName/$fileName", modes)
    catch
        run(`mkdir $dirName`)
        run(`touch $dirName/$fileName`)
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


function storeSamplePoints()
    xSampleFile = myOpen(;fileName="xSamples.csv", dirName=dirName, csvHeader="x(m)")
    ySampleFile = myOpen(;fileName="ySamples.csv", dirName=dirName, csvHeader="y(m)")
    zSampleFile = myOpen(;fileName="zSamples.csv", dirName=dirName, csvHeader="z(m)")
    for x in xs
        write(xSampleFile, "$(round(x, sigdigits=5))\n")
    end
    for y in ys
        write(ySampleFile, "$(round(y, sigdigits=5))\n")
    end
    for z in zs
        write(zSampleFile, "$(round(z, sigdigits=5))\n")
    end
    close(xSampleFile)
    close(ySampleFile)
    close(zSampleFile)
end


# Main

storeSamplePoints()

let resultsInZ0PlaneFile, bs, resultFile
	# init resultInZ0PlaneFile
    resultsInZ0PlaneFile = myOpen(;fileName="resultsInZ0Plane.csv", dirName=dirName, csvHeader=nothing)
	# init final results
    bs = zeros((length(xs), length(ys), length(zs)))
	# main calculation of b at sample points
    for (zIndex, zValue) in enumerate(zs), (yIndex, yValue) in enumerate(ys), (xIndex, xValue) in enumerate(xs)
        b::Float64 = @time calculateBAt(xValue, yValue, zValue)
        bs[xIndex, yIndex, zIndex] = b
        # result storation
        if zIndex == length(zs)
            if xIndex == length(xs)
                write(resultsInZ0PlaneFile, "$b\n")
            else
                write(resultsInZ0PlaneFile, "$b,")
            end
        end
    end
    close(resultsInZ0PlaneFile)

    meanB = mean(bs)
    minB = min(bs...)
    maxB = max(bs...)
    variationRate = (maxB - minB)/meanB
    println("variationRate = $variationRate")
    println("meanB = $meanB")

    resultFile = myOpen(;fileName="result.csv", dirName=dirName, csvHeader="variationRate(%),meanB(mT)")
    write(resultFile, "$(round(variationRate*100, sigdigits=6)),")
    write(resultFile, "$(round(meanB*1000, sigdigits=6))")
    close(resultFile)
end
