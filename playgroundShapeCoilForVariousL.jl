# Constants

const myu0 = 4pi*1e-7
const I = 100  # 1[A]
const N = 500
const l = 0.1  # 10cm
const h = 0.05  # 5cm
const R = h
const d = 0.58h


# Variables

const sourceIntervals = 100  # per part
const sampleIntervals = 500
const lower = 0.1  # 10cm
const upper = 1.0
const spacing = 0.01
const z = 0.1h

const sampleXHalfRanges = ( i*l for i=lower:spacing:upper )

const sampleYHalfRange = 0.1h
const sampleYSpacing = 2*sampleYHalfRange/sampleIntervals
const sampleYIntervals = [ -sampleYHalfRange + sampleYSpacing*i for i=1:sampleIntervals-1 ]


# Models

struct Point
    x::Float64
    y::Float64
    z::Float64
end
Point(array::Array{Float64, 1}) = Point(array[1], array[2], array[3])


function numericalIntegrateOf(f::Function; upperLimit::Float64, lowerLimit::Float64, n::Integer, point::Point)::Array{Float64, 1}
    sumOfIntervals = let
        intervals = [ lowerLimit + k*(upperLimit-lowerLimit)/n for k in 1:n-1 ]
        vector = [0.0; 0.0; 0.0]
        for k in intervals
            vector += f(;point=point, var=k)
        end
        vector
    end

    # intervals = [ lowerLimit + k*(upperLimit-lowerLimit)/n for k=1:n-1 ]
    # sumOfIntervals = sum(f.(;point=point, var=intervals) for k in intervals)

    return (upperLimit-lowerLimit)/n * ( (f(;point=point, var=upperLimit)+f(;point=point, var=lowerLimit))/2 + sumOfIntervals )
end


function c1(;point::Point, var::Float64)::Array{Float64, 1}
    xVector = 0
    yVector = -(point.z+d)
    zVector = (point.y+h)
    denominator = ( (point.x-var)^2 + (point.y+h)^2 + (point.z+d)^2 ) ^ (1.5)
    return [ xVector; yVector/denominator; zVector/denominator ]
end


function c2(;point::Point, var::Float64)::Array{Float64, 1}
    xVector = 0.0
    yVector = (point.z+d)
    zVector = -(point.y-h)
    denominator = ( (point.x-var)^2 + (point.y-h)^2 + (point.z+d)^2 ) ^ (1.5)
    return [ xVector; yVector/denominator; zVector/denominator ]
end


function c3(;point::Point, var::Float64)::Array{Float64, 1}
    xVector = (point.z+d)R*cos(var)
    yVector = (point.z+d)R*sin(var)
    zVector = -R*( (point.y-R*sin(var))sin(var) + (point.x-l-R*cos(var))cos(var) )
    denominator = ( (point.x-l-R*cos(var))^2 + (point.y-R*sin(var))^2 + (point.z+d)^2 ) ^ (1.5)
    return [ xVector/denominator; yVector/denominator; zVector/denominator ]
end


function c4(;point::Point, var::Float64)::Array{Float64, 1}
    xVector = (point.z+d)R*cos(var)
    yVector = (point.z+d)R*sin(var)
    zVector = -R*( (point.y-R*sin(var))sin(var) + (point.x+l-R*cos(var))cos(var) )
    denominator = ( (point.x+l-R*cos(var))^2 + (point.y-R*sin(var))^2 + (point.z+d)^2 ) ^ (1.5)
    return [ xVector/denominator; yVector/denominator; zVector/denominator ]
end


function c5(;point::Point, var::Float64)::Array{Float64, 1}
    xVector = 0
    yVector = -(point.z-d)
    zVector = (point.y+h)
    denominator = ( (point.x-var)^2 + (point.y+h)^2 + (point.z-d)^2 ) ^ (1.5)
    return [ xVector; yVector/denominator; zVector/denominator ]
end


function c6(;point::Point, var::Float64)::Array{Float64, 1}
    xVector = 0.0
    yVector = (point.z-d)
    zVector = -(point.y-h)
    denominator = ( (point.x-var)^2 + (point.y-h)^2 + (point.z-d)^2 ) ^ (1.5)
    return [ xVector; yVector/denominator; zVector/denominator ]
end


function c7(;point::Point, var::Float64)::Array{Float64, 1}
    xVector = (point.z-d)R*cos(var)
    yVector = (point.z-d)R*sin(var)
    zVector = -R*( (point.y-R*sin(var))sin(var) + (point.x-l-R*cos(var))cos(var) )
    denominator = ( (point.x-l-R*cos(var))^2 + (point.y-R*sin(var))^2 + (point.z-d)^2 ) ^ (1.5)
    return [ xVector/denominator; yVector/denominator; zVector/denominator ]
end


function c8(;point::Point, var::Float64)::Array{Float64, 1}
    xVector = (point.z-d)R*cos(var)
    yVector = (point.z-d)R*sin(var)
    zVector = -R*( (point.y-R*sin(var))sin(var) + (point.x+l-R*cos(var))cos(var) )
    denominator = ( (point.x+l-R*cos(var))^2 + (point.y-R*sin(var))^2 + (point.z-d)^2 ) ^ (1.5)
    return [ xVector/denominator; yVector/denominator; zVector/denominator ]
end


function BAtPointFromLower(point::Point)::Array{Float64, 1}
    result::Array{Float64, 1} = numericalIntegrateOf(c1; lowerLimit=-l, upperLimit=l, n=sourceIntervals, point=point)
    result += numericalIntegrateOf(c2; lowerLimit=-l, upperLimit=l, n=sourceIntervals, point=point)
    result += numericalIntegrateOf(c3; lowerLimit=-pi/2, upperLimit=pi/2, n=sourceIntervals, point=point)
    result += numericalIntegrateOf(c4; lowerLimit=pi/2, upperLimit=3pi/2, n=sourceIntervals, point=point)
    return result
end


function BAtPointFromUpper(point::Point)::Array{Float64, 1}
    result::Array{Float64, 1} = numericalIntegrateOf(c5; lowerLimit=-l, upperLimit=l, n=sourceIntervals, point=point)
    result += numericalIntegrateOf(c6; lowerLimit=-l, upperLimit=l, n=sourceIntervals, point=point)
    result += numericalIntegrateOf(c7; lowerLimit=-pi/2, upperLimit=pi/2, n=sourceIntervals, point=point)
    result += numericalIntegrateOf(c8; lowerLimit=pi/2, upperLimit=3pi/2, n=sourceIntervals, point=point)
    return result
end


using Statistics
function calculateBin(;planeZValue::Float64, sampleXIntervals::Array{Float64, 1})
    samplePoints = sampleIntervals-1
    result = map(zeros((samplePoints, samplePoints))) do x
        Point(x, x, x)
    end

    @time for (xIndex, xValue) in enumerate(sampleXIntervals), (yIndex, yValue) in enumerate(sampleYIntervals)
        resultInArray = myu0/(4pi)*(N*I) .* ( BAtPointFromLower(Point(xValue, yValue, planeZValue)) .+ BAtPointFromUpper(Point(xValue, yValue, planeZValue)) )
        result[xIndex, yIndex] = Point(resultInArray)
    end

    zValues = [ point.z for point in result ]
    minBOfZElement = min(zValues...)
    maxBOfZElement = max(zValues...)
    meanBOfZElement = mean(zValues)
    return (result, minBOfZElement, maxBOfZElement, meanBOfZElement)
end


function openWithNewDirIfNeeded(;dirName::String, fileName::String, modes::String)
    isSuccessful = true
    file = try
        open("$(dirName)/$fileName", modes)
    catch
        run(`mkdir $dirName`)
        run(`touch $dirName/$fileName`)
        isSuccessful = false
    end

    if isSuccessful == false
        file = open("$(dirName)/$fileName", modes)
    end
    return file
end


function storeVariationRatesOnly(file; index::Int, length::Float64, meanB::Float64, variationRate::Float64)
    # write variation rate
    if index == 1
        write(file, "length(h),variationRate,meanB[mT]\n")
    end
    write(file, "$(round(length/l, sigdigits=5)),$(round(variationRate, sigdigits=5)),$(round(meanB*1e3, sigdigits=4))\n")
end


# Main

const dirName = "I=$(round(I, sigdigits=2))_N=$(round(Int, N))_z=$(round(z/h, sigdigits=2))h_y=$(round(sampleYHalfRange/h, sigdigits=2))h_d=$(round(d/h, sigdigits=2))h"
const fileName = "variantionRateUnderVariousLFrom$(lower)To$(upper).csv"
const file = openWithNewDirIfNeeded(;dirName=dirName, fileName=fileName, modes="w")

for (index, sampleXHalfRange) in enumerate(sampleXHalfRanges)
    sampleXSpacing = 2*sampleXHalfRange/sampleIntervals
    sampleXIntervals = [ -sampleXHalfRange + sampleXSpacing*i for i=1:sampleIntervals-1 ]

    resultIn0Plane, minBOfZElementIn0Plane, maxBOfZElementIn0Plane, meanBOfZElementIn0Plane = calculateBin(planeZValue=0.0, sampleXIntervals=sampleXIntervals)
    resultInMiddlePlane, minBOfZElementInMiddlePlane, maxBOfZElementInMiddlePlane, meanBOfZElementInMiddlePlane = calculateBin(planeZValue=z/2, sampleXIntervals=sampleXIntervals)
    resultInUpperPlane, minBOfZElementInUpperPlane, maxBOfZElementInUpperPlane, meanBOfZElementInUpperPlane = calculateBin(planeZValue=z, sampleXIntervals=sampleXIntervals)

    minB = min(minBOfZElementIn0Plane, minBOfZElementInUpperPlane, minBOfZElementInMiddlePlane)
    maxB = max(maxBOfZElementIn0Plane, maxBOfZElementInUpperPlane, maxBOfZElementInMiddlePlane)
    meanB = mean((meanBOfZElementIn0Plane, meanBOfZElementInUpperPlane, meanBOfZElementInMiddlePlane))
    variationRate = (maxB-minB)/meanB*100

    println("Current Loop Area at 2h = $(2h*100)cm, 2d = $(round(2d*100, sigdigits=3))cm, 2l = $(2l*100)cm; ")
    println("Conducting Area at x = $(round(-sampleXHalfRange/l, sigdigits=4))l~$(round(sampleXHalfRange/l, sigdigits=2))l, y = $(round(sampleYHalfRange/h, sigdigits=2))h~$(round(sampleYHalfRange/h, sigdigits=2))h; samplePoint=$(sampleIntervals):" )
    println("min B of z elment when length = $(round(sampleXHalfRange/l, sigdigits=4))h is: $(minB*1e3) [mT]")
    println("max B of z elment when length = $(round(sampleXHalfRange/l, sigdigits=4))h is: $(maxB*1e3) [mT]")
    println("Magnetic Field Variance Rate: $variationRate%\n")

    global file
    storeVariationRatesOnly(file; index=index, length=sampleXHalfRange, meanB=meanB, variationRate=variationRate)
end

close(file)
