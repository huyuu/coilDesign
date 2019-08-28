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
@everywhere const h = 5e-2
# Measurement Area
@everywhere const X0 = 1e-2  # 1cm
@everywhere const Y0 = 1e-2  # 1cm
@everywhere const Z0 = 1e-2  # 1cm


# Variables

# Intervals into which a current source loop is cut, ex: c1. Should not be too small otherwise divergence condition could be raised."
@everywhere const sourceIntervals = 100
# Measurement points
@everywhere const sampleIntervals = 50
@everywhere const samplePoints = sampleIntervals+1
@everywhere const axisPoints = 100
# Coil Shape
@everywhere const dLowerCoeff = 0.3
@everywhere const dUpperCoeff = 0.9
@everywhere const lLowerCoeff = 1.0
@everywhere const lUpperCoeff = 5.0


# Children Variables

@everywhere const ds = LinRange(dLowerCoeff*h, dUpperCoeff*h, axisPoints)
@everywhere const dCoeffs = LinRange(dLowerCoeff, dUpperCoeff, axisPoints)
@everywhere const ls = LinRange(lLowerCoeff*h, lUpperCoeff*h, axisPoints)
@everywhere const lCoeffs = LinRange(lLowerCoeff, lUpperCoeff, axisPoints)
@everywhere const standardRadiusOfConductor = standardPhiOfConductor/2
@everywhere const conductorsPerLayer = div(thicknessOfGFRPWall, standardPhiOfConductor)
# sample points
@everywhere const xs = LinRange(-X0, X0, samplePoints)
@everywhere const ys = LinRange(-Y0, Y0, samplePoints)
@everywhere const zs = LinRange(-Z0, Z0, samplePoints)
# File Operation
const dirName = "I=$(round(I, sigdigits=2))_N=$(round(Int, N))_d=From$(dLowerCoeff)To$(dUpperCoeff)h_l=From$(lLowerCoeff)To$(lUpperCoeff)h_conductorPhi=$(round(standardPhiOfConductor*1e3, sigdigits=3))mm"


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
@everywhere ParsOfSource(n::Int, d::Float64) = let
    index::PositionIndex = PositionIndex(n)

    global h
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
    varRateVector = (maxBVector .- minBVector) ./ meanBVector
    ResultsOfMeanVarRate(meanBVector, varRateVector)
end


@everywhere function numericalIntegrateOf(f::Function; upperLimit::Float64, lowerLimit::Float64, n::Integer)::BVector
    sumOfIntervals = let
        temp::BVector = zeros(3)
        intervals = [ lowerLimit + k*(upperLimit-lowerLimit)/n for k in 1:n-1 ]
        for k in intervals
            temp .+ f(k)
        end
        temp
    end
    return (upperLimit-lowerLimit)/n * ( (f(upperLimit) .+ f(lowerLimit))/2 .+ sumOfIntervals )
end


@everywhere function c1(x::Float64, y::Float64, z::Float64; _x::Float64, _y::Float64, _z::Float64)::BVector
    xVector::Float64 = 0
    yVector::Float64 = -(z-_z)
    zVector::Float64 = y - _y
    denominator::Float64 = ( (x-_x)^2 + (y-_y)^2 + (z-_z)^2 ) ^ (1.5)
    return [xVector, yVector/denominator, zVector/denominator]
end


@everywhere function c2(x::Float64, y::Float64, z::Float64; _x::Float64, _y::Float64, _z::Float64)::BVector
    xVector::Float64 = 0
    yVector::Float64 = z-_z
    zVector::Float64 = -(y - _y)
    denominator::Float64 = ( (x-_x)^2 + (y-_y)^2 + (z-_z)^2 ) ^ (1.5)
    return [xVector, yVector/denominator, zVector/denominator]
end


@everywhere function c3(x::Float64, y::Float64, z::Float64; R::Float64, phi::Float64, _z::Float64, l::Float64)::BVector
    xVector::Float64 = R*cos(phi)*(z-_z)
    yVector::Float64 = R*sin(phi)*(z-_z)
    zVector::Float64 = -R*( (y-R*sin(phi))sin(phi) + (x-l-R*cos(phi))cos(phi) )
    denominator::Float64 = ( (x-l-R*cos(phi))^2 + (y-R*sin(phi))^2 + (z-_z)^2 ) ^ (1.5)
    return [xVector/denominator, yVector/denominator, zVector/denominator]
end


@everywhere function c4(x::Float64, y::Float64, z::Float64; R::Float64, phi::Float64, _z::Float64, l::Float64)::BVector
    xVector::Float64 = R*cos(phi)*(z-_z)
    yVector::Float64 = R*sin(phi)*(z-_z)
    zVector::Float64 = -R*( (y-R*sin(phi))sin(phi) + (x+l-R*cos(phi))cos(phi) )
    denominator::Float64 = ( (x+l-R*cos(phi))^2 + (y-R*sin(phi))^2 + (z-_z)^2 ) ^ (1.5)
    return [xVector/denominator, yVector/denominator, zVector/denominator]
end


@everywhere function BFromSingleTurn(x, y, z; n::Int, d::Float64, l::Float64)::BVector
    pars::ParsOfSource = ParsOfSource(n, d)
    local result::BVector = [0, 0, 0]
    # B from lower
    result .+= numericalIntegrateOf((_x) -> c1(x, y, z; _x=_x, _y=pars.c1_y, _z=pars.c1_z); lowerLimit=-l, upperLimit=l, n=sourceIntervals)
    result .+= numericalIntegrateOf((_x) -> c2(x, y, z; _x=_x, _y=pars.c1_y, _z=pars.c1_z); lowerLimit=-l, upperLimit=l, n=sourceIntervals)
    result .+= numericalIntegrateOf((phi) -> c3(x, y, z; R=pars.c3_R, phi=phi, _z=pars.c3_z, l=l); lowerLimit=-pi/2, upperLimit=pi/2, n=sourceIntervals)
    result .+= numericalIntegrateOf((phi) -> c4(x, y, z; R=pars.c4_R, phi=phi, _z=pars.c4_z, l=l); lowerLimit=pi/2, upperLimit=3pi/2, n=sourceIntervals)
    # B from upper
    result .+= numericalIntegrateOf((_x) -> c1(x, y, z; _x=_x, _y=pars.c5_y, _z=pars.c5_z); lowerLimit=-l, upperLimit=l, n=sourceIntervals)
    result .+= numericalIntegrateOf((_x) -> c2(x, y, z; _x=_x, _y=pars.c6_y, _z=pars.c6_z); lowerLimit=-l, upperLimit=l, n=sourceIntervals)
    result .+= numericalIntegrateOf((phi) -> c3(x, y, z; R=pars.c7_R, phi=phi, _z=pars.c7_z, l=l); lowerLimit=-pi/2, upperLimit=pi/2, n=sourceIntervals)
    result .+= numericalIntegrateOf((phi) -> c4(x, y, z; R=pars.c8_R, phi=phi, _z=pars.c8_z, l=l); lowerLimit=pi/2, upperLimit=3pi/2, n=sourceIntervals)
    return result
end


@everywhere function calculateBAtPoint(x::Float64, y::Float64, z::Float64; d::Float64, l::Float64)::BVector
    local results::BVector = zeros(3)

    map(1:N) do n
        results .+= BFromSingleTurn(x, y, z; n=n, d=d, l=l)
    end

    return μ/(4pi)*I * results
end


function calculateResultWhen(; d::Float64, l::Float64)::ResultsOfMeanVarRate
	# init results
    meanBVector = zeros(3)
    minBVector = ones(3)
    maxBVector = zeros(3)
	# for calculation
    totalSamplePoints = length(xs)*length(ys)*length(zs)
    futureBVectorsInCube::Array{Future} = []

    for (xIndex, xValue) in enumerate(xs), (yIndex, yValue) in enumerate(ys), (zIndex, zValue) in enumerate(zs)
        future = @spawn calculateBAtPoint(xValue, yValue, zValue; d=d, l=l)
        push!(futureBVectorsInCube, future)
    end

    map(futureBVectorsInCube) do futureBVector
        bVector = fetch(futureBVector)
        meanBVector .+= bVector
        map(1:3) do i
            minBVector[i] = bVector[i] < minBVector[i] ? bVector[i] : minBVector[i]
            maxBVector[i] = bVector[i] > maxBVector[i] ? bVector[i] : maxBVector[i]
        end
    end
    meanBVector ./= totalSamplePoints

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


function storeSamplePoints()
    dSampleFile = myOpen(;fileName="dSamples.csv", dirName=dirName, csvHeader="d(h)")
    lSampleFile = myOpen(;fileName="lSamples.csv", dirName=dirName, csvHeader="l(h)")
    for dCoeff in dCoeffs
        write(dSampleFile, "$(round(dCoeff, sigdigits=5))\n")
    end
    for lCoeff in lCoeffs
        write(lSampleFile, "$(round(lCoeff, sigdigits=5))\n")
    end
    close(dSampleFile)
    close(lSampleFile)
end


# Main

storeSamplePoints()

let
	for (dIndex, dValue) in enumerate(ds), (lIndex, lValue) in enumerate(ls)
        # calculation
        result = @time calculateResultWhen(; d=dValue, l=lValue)
		# open files
		meanBxFile = myOpen(;fileName="meanBx.csv", modes="a", dirName=dirName, csvHeader=nothing)
		meanByFile = myOpen(;fileName="meanBy.csv", modes="a", dirName=dirName, csvHeader=nothing)
		meanBzFile = myOpen(;fileName="meanBz.csv", modes="a", dirName=dirName, csvHeader=nothing)
		varRateXFile = myOpen(;fileName="varriationRateX.csv", modes="a", dirName=dirName, csvHeader=nothing)
		varRateYFile = myOpen(;fileName="varriationRateY.csv", modes="a", dirName=dirName, csvHeader=nothing)
		varRateZFile = myOpen(;fileName="varriationRateZ.csv", modes="a", dirName=dirName, csvHeader=nothing)
        # results storage
		if lIndex == length(ls)
			write(meanBxFile, "$(result.meanBVector[1])\n")
			write(meanByFile, "$(result.meanBVector[2])\n")
			write(meanBzFile, "$(result.meanBVector[3])\n")
            write(varRateXFile, "$(result.varRateVector[1])\n")
            write(varRateYFile, "$(result.varRateVector[2])\n")
            write(varRateZFile, "$(result.varRateVector[3])\n")
        else
            write(meanBxFile, "$(result.meanBVector[1]),")
            write(meanByFile, "$(result.meanBVector[2]),")
            write(meanBzFile, "$(result.meanBVector[3]),")
            write(varRateXFile, "$(result.varRateVector[1]),")
            write(varRateYFile, "$(result.varRateVector[2]),")
            write(varRateZFile, "$(result.varRateVector[3]),")
		end
        # closing files
        map((meanBxFile, meanByFile, meanBzFile, varRateXFile, varRateYFile, varRateZFile)) do file
            close(file)
        end
	end
end
