# This is the precise prediction for coil design.
using Distributed
addprocs(4)
@everywhere using Statistics


# Constants
# Physical Constants
@everywhere const μ = 4pi*1e-7
# Current
@everywhere const I = 3  # 1[A]
@everywhere const N = 200
# Ingredient
@everywhere const standardPhiOfConductor = 2.03e-3  # 2.03mm using AWG 18
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


# Variables

# Tolerable Variation Rate
@everywhere const standard = 0.01
# Error = 1%
@everywhere const ε = standard * 0.01
# Intervals into which a current source loop is cut, ex: c1. Should not be too small otherwise divergence condition could be raised."
@everywhere const sourceIntervals = 50
# Measurement points
@everywhere const sampleIntervals = 50
@everywhere const samplePoints = sampleIntervals+1
# File Operation
const dirName = "precise_I=$(round(I, sigdigits=2))_N=$(round(Int, N))_h=$(round(h*100, sigdigits=2))cm_X0=$(round(X0/h, sigdigits=2))h_Y0=$(round(Y0/h, sigdigits=2))h_Z0=$(round(Z0/h, sigdigits=2))h"


# Children Variables

@everywhere const xs = LinRange(-X0, X0, samplePoints)
@everywhere const ys = LinRange(-Y0, Y0, samplePoints)
@everywhere const zs = LinRange(-Z0, Z0, samplePoints)
@everywhere const standardRadiusOfConductor = standardPhiOfConductor/2
@everywhere const conductorsPerLayer = div(thicknessOfGFRPWall, standardPhiOfConductor)


# Models

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


@everywhere function numericalIntegrateOf(f::Function; upperLimit::Float64, lowerLimit::Float64, n::Integer)::Float64
    sumOfIntervals = let
        intervals = [ lowerLimit + k*(upperLimit-lowerLimit)/n for k in 1:n-1 ]
        sum( f(k) for k in intervals )
    end
    return (upperLimit-lowerLimit)/n * ( (f(upperLimit)+f(lowerLimit))/2 + sumOfIntervals )
end


@everywhere function c1(x::Float64, y::Float64, z::Float64; _x::Float64, _y::Float64, _z::Float64)::Float64
    zVector = y - _y
    denominator = ( (x-_x)^2 + (y-_y)^2 + (z-_z)^2 ) ^ (1.5)
    return zVector/denominator
end


@everywhere function c2(x::Float64, y::Float64, z::Float64; _x::Float64, _y::Float64, _z::Float64)::Float64
    zVector = -(y - _y)
    denominator = ( (x-_x)^2 + (y-_y)^2 + (z-_z)^2 ) ^ (1.5)
    return zVector/denominator
end


@everywhere function c3(x::Float64, y::Float64, z::Float64; R::Float64, phi::Float64, _z::Float64)::Float64
    zVector = -R*( (y-R*sin(phi))sin(phi) + (x-l-R*cos(phi))cos(phi) )
    denominator = ( (x-l-R*cos(phi))^2 + (y-R*sin(phi))^2 + (z-_z)^2 ) ^ (1.5)
    return zVector/denominator
end


@everywhere function c4(x::Float64, y::Float64, z::Float64; R::Float64, phi::Float64, _z::Float64)::Float64
    zVector = -R*( (y-R*sin(phi))sin(phi) + (x+l-R*cos(phi))cos(phi) )
    denominator = ( (x+l-R*cos(phi))^2 + (y-R*sin(phi))^2 + (z-_z)^2 ) ^ (1.5)
    return zVector/denominator
end


@everywhere function getSourceSegment(n::Int)::ParsOfSource
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


@everywhere function BFromSingleTurn(x, y, z; n::Int)
    pars::ParsOfSource = getSourceSegment(n)
    result::Float64 = 0
    # B from lower
    result += numericalIntegrateOf((_x) -> c1(x, y, z; _x=_x, _y=pars.c1_y, _z=pars.c1_z); lowerLimit=-l, upperLimit=l, n=sourceIntervals)
    result += numericalIntegrateOf((_x) -> c2(x, y, z; _x=_x, _y=pars.c1_y, _z=pars.c1_z); lowerLimit=-l, upperLimit=l, n=sourceIntervals)
    result += numericalIntegrateOf((phi) -> c3(x, y, z; R=pars.c3_R, phi=phi, _z=pars.c3_z); lowerLimit=-pi/2, upperLimit=pi/2, n=sourceIntervals)
    result += numericalIntegrateOf((phi) -> c4(x, y, z; R=pars.c4_R, phi=phi, _z=pars.c4_z); lowerLimit=pi/2, upperLimit=3pi/2, n=sourceIntervals)
    # B from upper
    result += numericalIntegrateOf((_x) -> c1(x, y, z; _x=_x, _y=pars.c5_y, _z=pars.c5_z); lowerLimit=-l, upperLimit=l, n=sourceIntervals)
    result += numericalIntegrateOf((_x) -> c2(x, y, z; _x=_x, _y=pars.c6_y, _z=pars.c6_z); lowerLimit=-l, upperLimit=l, n=sourceIntervals)
    result += numericalIntegrateOf((phi) -> c3(x, y, z; R=pars.c7_R, phi=phi, _z=pars.c7_z); lowerLimit=-pi/2, upperLimit=pi/2, n=sourceIntervals)
    result += numericalIntegrateOf((phi) -> c4(x, y, z; R=pars.c8_R, phi=phi, _z=pars.c8_z); lowerLimit=pi/2, upperLimit=3pi/2, n=sourceIntervals)
    return result
end


@everywhere function calculateBAt(x::Float64, y::Float64, z::Float64)::Float64
    local results::Float64 = 0
    futureResultsOfSingleTurn::Array{Future} = []

    map(1:N) do n
        future = @spawn BFromSingleTurn(x, y, z; n=n)
        push!(futureResultsOfSingleTurn, future)
    end

    map(futureResultsOfSingleTurn) do futureResult
        results += fetch(futureResult)
    end

    return results
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
    xSampleFile = myOpen(;fileName="xSamples", dirName=dirName, csvHeader="x(m)")
    ySampleFile = myOpen(;fileName="ySamples", dirName=dirName, csvHeader="y(m)")
    zSampleFile = myOpen(;fileName="zSamples", dirName=dirName, csvHeader="z(m)")
    for x in xs
        write(xSampleFile, "$(round(x, sigdigits=5))\n")
    end
    for y in ys
        write(ySampleFile, "$(round(y, sigdigits=5))\n")
    end
    for z in zs
        write(zSampleFile, "$(round(z, sigdigits=5))\n")
    end
end


# Main

storeSamplePoints()

resultsFile = myOpen(;fileName="resultsInZ0Plane", dirName=dirName, csvHeader=nothing)
results = zeros((length(xs), length(ys), length(zs)))
for (xIndex, xValue) in enumerate(xs), (yIndex, yValue) in enumerate(ys), (zIndex, zValue) in enumerate(zs)
    global results
    b::Float64 = @time calculateBAt(xValue, yValue, zValue)
    results[xIndex, yIndex, zIndex] = b
    if zIndex == length(zs)
        if xIndex != length(xs)
            write(resultsFile, "$b,")
        else
            write(resultsFile, "$b\n")
        end
    end
end
close(resultsFile)

let resultFile
    meanB = mean(results)
    minB = min(results...)
    maxB = max(results...)
    variationRate = (maxB - minB)/meanB

    resultFile = myOpen(;fileName="result", dirName=dirName, csvHeader="variationRate(%),meanB(mT)")
    write(resultFile, "$(round(variationRate*100, sigdigits=6))")
    write(resultFile, "$(round(meanB*1000, sigdigits=6))")
    close(resultFile)
end
