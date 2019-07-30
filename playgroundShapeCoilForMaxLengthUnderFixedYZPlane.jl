using Distributed

# Constants

@everywhere const myu0 = 4pi*1e-7
@everywhere const I = 50  # 1[A]
@everywhere const N = 50
@everywhere const h = 0.05  # 5cm
@everywhere const R = h
@everywhere const Y0 = 0.1h
@everywhere const Z0 = 0.1h


# Variables

# Intervals into which a line is cut, ex: c1. Should not be too small otherwise divergence condition could be raised."
@everywhere const sourceIntervals = 50
@everywhere const sampleIntervals = 50
@everywhere const samplePoints = sampleIntervals+1
# Points for the outer most loop, therefore should not be too much."
@everywhere const axisPoints = 50


@everywhere const ds = let
    lowerCoeff = 0.1
    upperCoeff = 1.0
    n::Int = axisPoints
    ( i*h for i=range(lowerCoeff, stop=upperCoeff, length=n) )
end

@everywhere const ls = let
    lowerCoeff = 1.0
    upperCoeff = 5.0
    n::Int = axisPoints
    ( i*h for i=range(lowerCoeff, stop=upperCoeff, length=n) )
end

# const sampleXHalfRanges = ( i*l for i=lower:spacing:upper )
#
# const sampleYHalfRange = 0.1h
# const sampleYSpacing = 2*sampleYHalfRange/sampleIntervals
# const sampleYIntervals = [ -sampleYHalfRange + sampleYSpacing*i for i=1:sampleIntervals-1 ]


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


using Statistics
function calculateVariationRateAndMeanBWhen(X0; d::Float64, l::Float64)::Tuple{Float64, Float64}
    result = zeros((samplePoints, samplePoints, samplePoints))
    # Set sample points in Y-Z plane
    ySamplePoints = LinRange(-Y0, Y0, samplePoints)
    zSamplePoints = LinRange(-Z0, Z0, samplePoints)
    xSamplePoints = LinRange(0.0, X0, round(Int, sampleIntervals/2 + 1))
    allPoints = length(xSamplePoints) * length(ySamplePoints) * length(zSamplePoints)

    for (xIndex, xValue) in enumerate(xSamplePoints), (yIndex, yValue) in enumerate(ySamplePoints), (zIndex, zValue) in enumerate(zSamplePoints)
        resultFromUpper = BAtPointFromUpper(Point(xValue, yValue, zValue); d=d, l=l)
        resultFromLower = BAtPointFromLower(Point(xValue, yValue, zValue); d=d, l=l)
        result[xIndex, yIndex, zIndex] = myu0/(4pi)*(N*I) .* ( resultFromLower .+ resultFromUpper )
        # resultInArray::Array{Float64, 1} = myu0/(4pi)*(N*I) .* ( BAtPointFromLower(Point(xValue, yValue, zValue); d=d, l=l) .+ BAtPointFromUpper(Point(xValue, yValue, zValue); d=d, l=l) )
    end
    meanBOfZElement = mean(result)
    minBOfZElement = min(result...)
    maxBOfZElement = max(result...)

    variationRate = (maxBOfZElement-minBOfZElement)/meanBOfZElement
    return variationRate, meanBOfZElement
end


function solveByInterpolation(;xLower::Float64, xUpper::Float64, dValue::Float64, lValue::Float64)::Tuple{Float64, Float64, Float64}
    standard = 0.0
    signAndResultAt(x) = let
        variationRate, meanB = calculateVariationRateAndMeanBWhen(x; d=dValue, l=lValue)
        signOfResult = sign(variationRate-standard)
        (signOfResult, variationRate, meanB)
    end
    xLowerSign, xLowerVarRate, xLowerMeanB = signAndResultAt(xLower)
    xUpperSign, xUpperVarRate, xUpperMeanB = signAndResultAt(xUpper)
    xMiddleSign, xMiddleVarRate, xMiddleMeanB = signAndResultAt((xUpper+xLower)/2)

    while true
        if xLowerSign == xUpperSign == xMiddleSign
            return (standard, 1.0, 0.0)  # no solution
        elseif isCloseEnough(xUpperSign, standard)
            return (xUpper, xUpperVarRate, xUpperMeanB)
        elseif isCloseEnough(xMiddleSign, standard)
            return (xMiddle, xMiddleVarRate, xMiddleMeanB)
        elseif isCloseEnough(xLowerSign, standard)
            return (xLower, xLowerVarRate, xLowerMeanB)

        elseif xLowerSign == xMiddleSign != xUpperSign
            xLower = xMiddle
            xLowerSign, xLowerVarRate, xLowerMeanB = signAndResultAt(xLower)
            xMiddle = (xMiddle+xUpper)/2
            xMiddleSign, xMiddleVarRate, xMiddleMeanB = signAndResultAt(xMiddle)
            continue
        elseif xLowerSign != xMiddleSign == xUpperSign
            xUpper = xMiddle
            xUpperSign, xUpperVarRate, xUpperMeanB = signAndResultAt(xUpper)
            xMiddle = (xLower+xMiddle)/2
            xMiddleSign, xMiddleVarRate, xMiddleMeanB = signAndResultAt(xMiddle)
            continue
        else
            error("Here")
        end
    end
end


function isCloseEnough(a::Number, b::Float64; ε::Float64=ε)
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
        write(file, csvHeader)
    else  # csvHeader == nothing
        error("csvHeader shouldn't be nothing.")
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
    fileX0s = myOpen(;fileName="maxX0s_D=$(dsLower)To$(dsUpper)_L=$(lsLower)To$(lsUpper).csv", dirName=dirName, withNewDirIfNeeded=true, csvHeader="maxX0(h)")
    fileMeanBs = myOpen(;fileName="meanBs_D=$(dsLower)To$(dsUpper)_L=$(lsLower)To$(lsUpper).csv", dirName=dirName, withNewDirIfNeeded=true, csvHeader="meanBs[mT]")

    return Dict("d"=>fileD, "l"=>fileL, "X0s"=>fileX0s, "meanBs"=>fileMeanBs)
end


function closeFiles(files::Dict{String, IOStream})
    map(values(files)) do file
        close(file)
    end
end


function storeX0AndMeanBs(files::Dict{String, IOStream}; shouldEndLine::Bool, X0::Float64, meanB::Float64)
    if shouldEndLine
        write(files["X0s"], round(X0, sigdigits=5), "\n")
        write(files["meanBs"], round(meanB, sigdigits=5), "\n")
    else
        write(files["X0s"], round(X0, sigdigits=5), ",")
        write(files["meanBs"], round(meanB, sigdigits=5), ",")
    end
end


function storeSamplePoints(files::Dict{String, IOStream})
    map(ds, ls) do d, l
        write(files["d"], round(d/h, sigdigits=4), "\n")
        write(files["l"], round(l/h, sigdigits=4), "\n")
    end
end


# Main
const files = openFiles()
const ε = 0.01 * 0.1  # error = 10%

for (dIndex, dValue) in enumerate(ds), (lIndex, lValue) in enumerate(ls)
    X0::Float64, variationRate::Float64, meanB::Float64 = @time solveByInterpolation(; xLower=0.0, xUpper=copy(lValue), dValue=dValue, lValue=lValue)
    println("(d:$(round(dValue/h, sigdigits=4))h, l:$(round(lValue/h, sigdigits=4))h): maxX0 = $(round(X0/h, sigdigits=6))h")
    storeX0AndMeanBs(files; shouldEndLine=lIndex==length(ls), X0=X0, meanB=meanB)
end
storeSamplePoints(files)




closeFiles(files)
