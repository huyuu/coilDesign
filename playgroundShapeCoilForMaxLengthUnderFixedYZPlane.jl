using Distributed
addprocs(6)
@everywhere using Statistics


# Constants

@everywhere const myu0 = 4pi*1e-7
@everywhere const I = 100  # 1[A]
@everywhere const N = 400
@everywhere const h = 0.05  # 5cm
@everywhere const R = h
@everywhere const Y0 = 0.1h
@everywhere const Z0 = 0.1h


# Variables

@everywhere const standard = 0.01  # tolerable variation rate
@everywhere const ε = standard * 0.01  # error = 1%
# Intervals into which a line is cut, ex: c1. Should not be too small otherwise divergence condition could be raised."
@everywhere const sourceIntervals = 50
@everywhere const sampleIntervals = 50
@everywhere const samplePoints = sampleIntervals+1
# Points for the outer most loop, therefore should not be too much."
@everywhere const axisPoints = 300


@everywhere const ds = let
    lowerCoeff = 0.3
    upperCoeff = 0.9
    n::Int = axisPoints
    ( i*h for i=range(lowerCoeff, stop=upperCoeff, length=n) )
end

@everywhere const ls = let
    lowerCoeff = 2.0
    upperCoeff = 10.0
    n::Int = axisPoints
    ( i*h for i=range(lowerCoeff, stop=upperCoeff, length=n) )
end


# Models

@everywhere struct Point
    x::Float64
    y::Float64
    z::Float64
end
@everywhere Point(array::Array{Float64, 1}) = Point(array[1], array[2], array[3])


@everywhere function numericalIntegrateOf(f::Function; upperLimit::Float64, lowerLimit::Float64, d::Float64, l::Float64, n::Integer, point::Point)::Float64
    sumOfIntervals = let
        intervals = [ lowerLimit + k*(upperLimit-lowerLimit)/n for k in 1:n-1 ]
        sum( f(;point=point, var=k, d=d, l=l) for k in intervals )
    end

    # intervals = [ lowerLimit + k*(upperLimit-lowerLimit)/n for k=1:n-1 ]
    # sumOfIntervals = sum(f.(;point=point, var=intervals) for k in intervals)

    return (upperLimit-lowerLimit)/n * ( (f(;point=point, var=upperLimit, d=d, l=l)+f(;point=point, var=lowerLimit, d=d, l=l))/2 + sumOfIntervals )
end


@everywhere function c1(;point::Point, var::Float64, d::Float64, l::Float64)::Float64
    xVector = 0
    yVector = -(point.z+d)
    zVector = (point.y+h)
    denominator = ( (point.x-var)^2 + (point.y+h)^2 + (point.z+d)^2 ) ^ (1.5)
    # return [ xVector; yVector/denominator; zVector/denominator ]
    return zVector/denominator
end


@everywhere function c2(;point::Point, var::Float64, d::Float64, l::Float64)::Float64
    xVector = 0.0
    yVector = (point.z+d)
    zVector = -(point.y-h)
    denominator = ( (point.x-var)^2 + (point.y-h)^2 + (point.z+d)^2 ) ^ (1.5)
    # return [ xVector; yVector/denominator; zVector/denominator ]
    return zVector/denominator
end


@everywhere function c3(;point::Point, var::Float64, d::Float64, l::Float64)::Float64
    xVector = (point.z+d)R*cos(var)
    yVector = (point.z+d)R*sin(var)
    zVector = -R*( (point.y-R*sin(var))sin(var) + (point.x-l-R*cos(var))cos(var) )
    denominator = ( (point.x-l-R*cos(var))^2 + (point.y-R*sin(var))^2 + (point.z+d)^2 ) ^ (1.5)
    # return [ xVector/denominator; yVector/denominator; zVector/denominator ]
    return zVector/denominator
end


@everywhere function c4(;point::Point, var::Float64, d::Float64, l::Float64)::Float64
    xVector = (point.z+d)R*cos(var)
    yVector = (point.z+d)R*sin(var)
    zVector = -R*( (point.y-R*sin(var))sin(var) + (point.x+l-R*cos(var))cos(var) )
    denominator = ( (point.x+l-R*cos(var))^2 + (point.y-R*sin(var))^2 + (point.z+d)^2 ) ^ (1.5)
    # return [ xVector/denominator; yVector/denominator; zVector/denominator ]
    return zVector/denominator
end


@everywhere function c5(;point::Point, var::Float64, d::Float64, l::Float64)::Float64
    xVector = 0
    yVector = -(point.z-d)
    zVector = (point.y+h)
    denominator = ( (point.x-var)^2 + (point.y+h)^2 + (point.z-d)^2 ) ^ (1.5)
    # return [ xVector; yVector/denominator; zVector/denominator ]
    return zVector/denominator
end


@everywhere function c6(;point::Point, var::Float64, d::Float64, l::Float64)::Float64
    xVector = 0.0
    yVector = (point.z-d)
    zVector = -(point.y-h)
    denominator = ( (point.x-var)^2 + (point.y-h)^2 + (point.z-d)^2 ) ^ (1.5)
    # return [ xVector/denominator; yVector/denominator; zVector/denominator ]
    return zVector/denominator
end


@everywhere function c7(;point::Point, var::Float64, d::Float64, l::Float64)::Float64
    xVector = (point.z-d)R*cos(var)
    yVector = (point.z-d)R*sin(var)
    zVector = -R*( (point.y-R*sin(var))sin(var) + (point.x-l-R*cos(var))cos(var) )
    denominator = ( (point.x-l-R*cos(var))^2 + (point.y-R*sin(var))^2 + (point.z-d)^2 ) ^ (1.5)
    # return [ xVector/denominator; yVector/denominator; zVector/denominator ]
    return zVector/denominator
end


@everywhere function c8(;point::Point, var::Float64, d::Float64, l::Float64)::Float64
    xVector = (point.z-d)R*cos(var)
    yVector = (point.z-d)R*sin(var)
    zVector = -R*( (point.y-R*sin(var))sin(var) + (point.x+l-R*cos(var))cos(var) )
    denominator = ( (point.x+l-R*cos(var))^2 + (point.y-R*sin(var))^2 + (point.z-d)^2 ) ^ (1.5)
    # return [ xVector/denominator; yVector/denominator; zVector/denominator ]
    return zVector/denominator
end


@everywhere function BAtPointFromLower(point::Point; d::Float64, l::Float64)::Float64
    result::Float64 = numericalIntegrateOf(c1; lowerLimit=-l, upperLimit=l, d=d, l=l, n=sourceIntervals, point=point)
    result += numericalIntegrateOf(c2; lowerLimit=-l, upperLimit=l, d=d, l=l, n=sourceIntervals, point=point)
    result += numericalIntegrateOf(c3; lowerLimit=-pi/2, upperLimit=pi/2, d=d, l=l, n=sourceIntervals, point=point)
    result += numericalIntegrateOf(c4; lowerLimit=pi/2, upperLimit=3pi/2, d=d, l=l, n=sourceIntervals, point=point)
    return result
end


@everywhere function BAtPointFromUpper(point::Point; d::Float64, l::Float64)::Float64
    result::Float64 = numericalIntegrateOf(c5; lowerLimit=-l, upperLimit=l, d=d, l=l, n=sourceIntervals, point=point)
    result += numericalIntegrateOf(c6; lowerLimit=-l, upperLimit=l, d=d, l=l, n=sourceIntervals, point=point)
    result += numericalIntegrateOf(c7; lowerLimit=-pi/2, upperLimit=pi/2, d=d, l=l, n=sourceIntervals, point=point)
    result += numericalIntegrateOf(c8; lowerLimit=pi/2, upperLimit=3pi/2, d=d, l=l, n=sourceIntervals, point=point)
    return result
end


@everywhere function calculateVariationRateAndMeanBWhen(X0; d::Float64, l::Float64)::Tuple{Float64, Float64}
    # Set sample points in Y-Z plane
    ySamplePoints = LinRange(-Y0, Y0, samplePoints)
    zSamplePoints = LinRange(-Z0, Z0, samplePoints)
    xSamplePoints = LinRange(0.0, X0, round(Int, samplePoints/2))
    results = zeros((length(xSamplePoints), length(ySamplePoints), length(zSamplePoints)))

    for (xIndex, xValue) in enumerate(xSamplePoints), (yIndex, yValue) in enumerate(ySamplePoints), (zIndex, zValue) in enumerate(zSamplePoints)
        resultsFromUpper =  BAtPointFromUpper(Point(xValue, yValue, zValue); d=d, l=l)
        resultsFromLower = BAtPointFromLower(Point(xValue, yValue, zValue); d=d, l=l)
        results[xIndex, yIndex, zIndex] = myu0/(4pi)*(N*I) * ( resultsFromLower + resultsFromUpper )
        # resultInArray::Array{Float64, 1} = myu0/(4pi)*(N*I) .* ( BAtPointFromLower(Point(xValue, yValue, zValue); d=d, l=l) .+ BAtPointFromUpper(Point(xValue, yValue, zValue); d=d, l=l) )
    end
    # println(errorDescription)
    meanBOfZElement = mean(results)
    minBOfZElement = min(results...)
    maxBOfZElement = max(results...)

    variationRate = (maxBOfZElement-minBOfZElement)/meanBOfZElement
    # println("X0 = $X0: min = $minBOfZElement, max = $maxBOfZElement, meanB = $meanBOfZElement, var = $variationRate")
    return (variationRate, meanBOfZElement)
end


function solveByInterpolation(;xLowerOrigin::Float64, xUpperOrigin::Float64, dValue::Float64, lValue::Float64)::Tuple{Float64, Float64, Float64}
    signAndResultAt(x) = let
        variationRate, meanB = calculateVariationRateAndMeanBWhen(x; d=dValue, l=lValue)
        signOfResult = sign(variationRate-standard)
        (signOfResult, variationRate, meanB)
    end
    xLower::Float64 = copy(xLowerOrigin)
    xUpper::Float64 = copy(xUpperOrigin)
    xMiddle = (xUpper+xLower)/2

    xLowerFutureResults = @spawn signAndResultAt(xLower)
    xUpperFutureResults = @spawn signAndResultAt(xUpper)
    xMiddleFutureResults = @spawn signAndResultAt(xMiddle)

    xLowerSign::Int, xLowerVarRate::Float64, xLowerMeanB::Float64 = fetch(xLowerFutureResults)
    xUpperSign::Int, xUpperVarRate::Float64, xUpperMeanB::Float64 = fetch(xUpperFutureResults)
    xMiddleSign::Int, xMiddleVarRate::Float64, xMiddleMeanB::Float64 = fetch(xMiddleFutureResults)

    while true
        if xLowerSign == xUpperSign == xMiddleSign
            println("return at first calculate with xLowerVarRate: $(xLowerVarRate), xMiddleVarRate: $(xMiddleVarRate), xUpperVarRate: $(xUpperVarRate)")
            return (0.0, 1.0, 0.0)  # no solution
        elseif isCloseEnough(xUpperVarRate, standard)
            return (xUpper, xUpperVarRate, xUpperMeanB)
        elseif isCloseEnough(xMiddleVarRate, standard)
            return (xMiddle, xMiddleVarRate, xMiddleMeanB)
        elseif isCloseEnough(xLowerVarRate, standard)
            return (xLower, xLowerVarRate, xLowerMeanB)

        elseif xLowerSign == xMiddleSign != xUpperSign
            xLower = xMiddle
            xLowerFutureResults = @spawn signAndResultAt(xLower)

            xMiddle = (xMiddle+xUpper)/2
            xMiddleFutureResults = @spawn signAndResultAt(xMiddle)

            xLowerSign, xLowerVarRate, xLowerMeanB = fetch(xLowerFutureResults)
            xMiddleSign, xMiddleVarRate, xMiddleMeanB = fetch(xMiddleFutureResults)

        elseif xLowerSign != xMiddleSign == xUpperSign
            xUpper = xMiddle
            xUpperFutureResults = @spawn signAndResultAt(xUpper)

            xMiddle = (xLower+xMiddle)/2
            xMiddleFutureResults = @spawn signAndResultAt(xMiddle)

            xUpperSign, xUpperVarRate, xUpperMeanB = fetch(xUpperFutureResults)
            xMiddleSign, xMiddleVarRate, xMiddleMeanB = fetch(xMiddleFutureResults)

        else
            println("$(xLowerSign), $xMiddleSign, $xUpperSign")
            error("Here")
        end
    end
end


function isCloseEnough(a::Number, b::Number; ε::Float64=ε)::Bool
    abs(a - b) < ε
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


function openFiles()::Dict{String, IOStream}
    dirName = "I=$(round(I, sigdigits=2))_N=$(round(Int, N))_h=$(round(h*100, sigdigits=2))cm_Y0=Z0=$(round(Y0/h, sigdigits=2))h"
    dsArray = collect(ds)
    lsArray = collect(ls)
    dsLower = "$(round(dsArray[1]/h, sigdigits=2))h"
    dsUpper = "$(round(dsArray[end]/h, sigdigits=2))h"
    lsLower = "$(round(lsArray[1]/h, sigdigits=2))h"
    lsUpper = "$(round(lsArray[end]/h, sigdigits=2))h"

    fileD = myOpen(;fileName="dSamplesFrom$(dsLower)To$(dsUpper).csv", dirName=dirName, withNewDirIfNeeded=true, csvHeader="d(h)")
    fileL = myOpen(;fileName="lSamplesFrom$(lsLower)To$(lsUpper).csv", dirName=dirName, withNewDirIfNeeded=true, csvHeader="l(h)")
    fileX0s = myOpen(;fileName="maxX0s_D=$(dsLower)To$(dsUpper)_L=$(lsLower)To$(lsUpper).csv", dirName=dirName, withNewDirIfNeeded=true, csvHeader=nothing)
    fileMeanBs = myOpen(;fileName="meanBs_D=$(dsLower)To$(dsUpper)_L=$(lsLower)To$(lsUpper).csv", dirName=dirName, withNewDirIfNeeded=true, csvHeader=nothing)

    return Dict("d"=>fileD, "l"=>fileL, "X0s"=>fileX0s, "meanBs"=>fileMeanBs)
end


function closeFiles(files::Dict{String, IOStream})
    map(values(files)) do file
        close(file)
    end
end


function storeX0AndMeanBs(files::Dict{String, IOStream}; shouldEndLine::Bool, X0::Float64, meanB::Float64)
    if shouldEndLine
        write(files["X0s"], "$(round(X0, sigdigits=5))\n")
        write(files["meanBs"], "$(round(meanB, sigdigits=5))\n")
    else
        write(files["X0s"], "$(round(X0, sigdigits=5)),")
        write(files["meanBs"], "$(round(meanB, sigdigits=5)),")
    end
end


function storeSamplePoints(files::Dict{String, IOStream})
    map(ds, ls) do d, l
        write(files["d"], "$(round(d/h, sigdigits=4))\n")
        write(files["l"], "$(round(l/h, sigdigits=4))\n")
    end
end


# Main
const files = openFiles()

storeSamplePoints(files)
for (dIndex, dValue) in enumerate(ds), (lIndex, lValue) in enumerate(ls)
    X0::Float64, variationRate::Float64, meanB::Float64 = @time solveByInterpolation(; xLowerOrigin=0.01ε, xUpperOrigin=lValue, dValue=dValue, lValue=lValue)
    println("(d:$(round(dValue/h, sigdigits=4))h, l:$(round(lValue/h, sigdigits=4))h): maxX0 = $(round(X0/h, sigdigits=6))h")
    storeX0AndMeanBs(files; shouldEndLine=lIndex==length(ls), X0=X0, meanB=meanB)
end





closeFiles(files)
