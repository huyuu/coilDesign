# Constants

const myu0 = 4pi*1e-7
const I = 100  # 1[A]
const N = 500
const h = 0.05  # 5cm
const l = 0.10  # 10cm
const R = h


# Variables

const sourceIntervals = 500  # per part
const sampleIntervals = 500
const ds = ( 0.01i*h for i=1:100 )
const z = 0.1h

const sampleXHalfRange = 0.1l
const sampleXSpacing = 2*sampleXHalfRange/sampleIntervals
const sampleXIntervals = [ -sampleXHalfRange + sampleXSpacing*i for i=1:sampleIntervals-1 ]

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


function numericalIntegrateOf(f::Function; upperLimit::Float64, lowerLimit::Float64, n::Integer, point::Point, distance::Float64)::Array{Float64, 1}
    sumOfIntervals = let
        intervals = [ lowerLimit + k*(upperLimit-lowerLimit)/n for k in 1:n-1 ]
        vector = [0.0; 0.0; 0.0]
        for k in intervals
            vector += f(;point=point, var=k, d=distance)
        end
        vector
    end

    # intervals = [ lowerLimit + k*(upperLimit-lowerLimit)/n for k=1:n-1 ]
    # sumOfIntervals = sum(f.(;point=point, var=intervals) for k in intervals)

    return (upperLimit-lowerLimit)/n * ( (f(;point=point, var=upperLimit, d=distance)+f(;point=point, var=lowerLimit, d=distance))/2 + sumOfIntervals )
end


function c1(;point::Point, var::Float64, d::Float64)::Array{Float64, 1}
    xVector = 0
    yVector = -(point.z+d)
    zVector = (point.y+h)
    denominator = ( (point.x-var)^2 + (point.y+h)^2 + (point.z+d)^2 ) ^ (1.5)
    return [ xVector; yVector/denominator; zVector/denominator ]
end


function c2(;point::Point, var::Float64, d::Float64)::Array{Float64, 1}
    xVector = 0.0
    yVector = (point.z+d)
    zVector = -(point.y-h)
    denominator = ( (point.x-var)^2 + (point.y-h)^2 + (point.z+d)^2 ) ^ (1.5)
    return [ xVector; yVector/denominator; zVector/denominator ]
end


function c3(;point::Point, var::Float64, d::Float64)::Array{Float64, 1}
    xVector = (point.z+d)R*cos(var)
    yVector = (point.z+d)R*sin(var)
    zVector = -R*( (point.y-R*sin(var))sin(var) + (point.x-l-R*cos(var))cos(var) )
    denominator = ( (point.x-l-R*cos(var))^2 + (point.y-R*sin(var))^2 + (point.z+d)^2 ) ^ (1.5)
    return [ xVector/denominator; yVector/denominator; zVector/denominator ]
end


function c4(;point::Point, var::Float64, d::Float64)::Array{Float64, 1}
    xVector = (point.z+d)R*cos(var)
    yVector = (point.z+d)R*sin(var)
    zVector = -R*( (point.y-R*sin(var))sin(var) + (point.x+l-R*cos(var))cos(var) )
    denominator = ( (point.x+l-R*cos(var))^2 + (point.y-R*sin(var))^2 + (point.z+d)^2 ) ^ (1.5)
    return [ xVector/denominator; yVector/denominator; zVector/denominator ]
end


function c5(;point::Point, var::Float64, d::Float64)::Array{Float64, 1}
    xVector = 0
    yVector = -(point.z-d)
    zVector = (point.y+h)
    denominator = ( (point.x-var)^2 + (point.y+h)^2 + (point.z-d)^2 ) ^ (1.5)
    return [ xVector; yVector/denominator; zVector/denominator ]
end


function c6(;point::Point, var::Float64, d::Float64)::Array{Float64, 1}
    xVector = 0.0
    yVector = (point.z-d)
    zVector = -(point.y-h)
    denominator = ( (point.x-var)^2 + (point.y-h)^2 + (point.z-d)^2 ) ^ (1.5)
    return [ xVector; yVector/denominator; zVector/denominator ]
end


function c7(;point::Point, var::Float64, d::Float64)::Array{Float64, 1}
    xVector = (point.z-d)R*cos(var)
    yVector = (point.z-d)R*sin(var)
    zVector = -R*( (point.y-R*sin(var))sin(var) + (point.x-l-R*cos(var))cos(var) )
    denominator = ( (point.x-l-R*cos(var))^2 + (point.y-R*sin(var))^2 + (point.z-d)^2 ) ^ (1.5)
    return [ xVector/denominator; yVector/denominator; zVector/denominator ]
end


function c8(;point::Point, var::Float64, d::Float64)::Array{Float64, 1}
    xVector = (point.z-d)R*cos(var)
    yVector = (point.z-d)R*sin(var)
    zVector = -R*( (point.y-R*sin(var))sin(var) + (point.x+l-R*cos(var))cos(var) )
    denominator = ( (point.x+l-R*cos(var))^2 + (point.y-R*sin(var))^2 + (point.z-d)^2 ) ^ (1.5)
    return [ xVector/denominator; yVector/denominator; zVector/denominator ]
end


function BAtPointFromLower(point::Point; distance::Float64)::Array{Float64, 1}
    result::Array{Float64, 1} = numericalIntegrateOf(c1; lowerLimit=-l, upperLimit=l, n=sourceIntervals, point=point, distance=distance)
    result += numericalIntegrateOf(c2; lowerLimit=-l, upperLimit=l, n=sourceIntervals, point=point, distance=distance)
    result += numericalIntegrateOf(c3; lowerLimit=-pi/2, upperLimit=pi/2, n=sourceIntervals, point=point, distance=distance)
    result += numericalIntegrateOf(c4; lowerLimit=pi/2, upperLimit=3pi/2, n=sourceIntervals, point=point, distance=distance)
    return result
end


function BAtPointFromUpper(point::Point; distance::Float64)::Array{Float64, 1}
    result::Array{Float64, 1} = numericalIntegrateOf(c5; lowerLimit=-l, upperLimit=l, n=sourceIntervals, point=point, distance=distance)
    result += numericalIntegrateOf(c6; lowerLimit=-l, upperLimit=l, n=sourceIntervals, point=point, distance=distance)
    result += numericalIntegrateOf(c7; lowerLimit=-pi/2, upperLimit=pi/2, n=sourceIntervals, point=point, distance=distance)
    result += numericalIntegrateOf(c8; lowerLimit=pi/2, upperLimit=3pi/2, n=sourceIntervals, point=point, distance=distance)
    return result
end


using Statistics
function calculateBin(;planeZValue::Float64, distance::Float64)
    samplePoints = sampleIntervals-1
    result = map(zeros((samplePoints, samplePoints))) do x
        Point(x, x, x)
    end

    @time for (xIndex, xValue) in enumerate(sampleXIntervals), (yIndex, yValue) in enumerate(sampleYIntervals)
        resultInArray = myu0/(4pi)*(N*I) .* ( BAtPointFromLower(Point(xValue, yValue, planeZValue); distance=distance) .+ BAtPointFromUpper(Point(xValue, yValue, planeZValue); distance=distance) )
        result[xIndex, yIndex] = Point(resultInArray)
    end

    zValues = [ point.z for point in result ]
    minBOfZElement = min(zValues...)
    maxBOfZElement = max(zValues...)
    meanBOfZElement = mean(zValues)
    return (result, minBOfZElement, maxBOfZElement, meanBOfZElement)
end


function store(result::Array{Point, 2}; index::Int, planeZValue::Float64, distance::Float64, maxB::Float64, minB::Float64, meanB::Float64, variationRate::Float64)
    # open files
    fileX, fileY, fileZ = map(["x", "y", "z"]) do var
        try
            open("$(dirName)/$(var)ElementsOfBWhenD=$(round(distance/h, sigdigits=3))h.csv", "w")
        catch errorStatement
            run(`mkdir $dirName`)
        finally
            open("$(dirName)/$(var)ElementsOfBWhenD=$(round(distance/h, sigdigits=3))h.csv", "w")
        end
    end
    fileXPoints, fileYPoints = map(('x', 'y')) do var
        open("$(dirName)/$(var)SamplePointsWhenD=$(round(distance/h, sigdigits=3))h.csv", "w")
    end
    fileVariationRate = try
        open("$(dirName)/variantionRateUnderVariousD.csv", "a")
    catch
        run(`touch $dirName/variantionRateUnderVariousD.csv`)
    finally
        open("$(dirName)/variantionRateUnderVariousD.csv", "a")
    end

    # write results to files
    for (xIndex, xValue) in enumerate(sampleXIntervals)
        for (yIndex, yValue) in enumerate(sampleYIntervals)
            point = result[xIndex, yIndex]
            map( zip((fileX, fileY, fileZ), (point.x, point.y, point.z)) ) do (file, value)
                write(file, "$(round(value, sigdigits=12)),")
            end
        end
        map((fileX, fileY, fileZ)) do file
            write(file, "\n")
        end
    end

    # write sample points
    map(sampleXIntervals) do value
        write(fileXPoints, "$(round(value, sigdigits=4))\n")
    end
    map(sampleYIntervals) do value
        write(fileYPoints, "$(round(value, sigdigits=4))\n")
    end
    # write variation rate
    index == 1 ? write(fileVariationRate, "d(h),variationRate,meanB[mT]\n") : nothing
    write(fileVariationRate, "$(round(distance/h, sigdigits=4)),$(round(variationRate, sigdigits=5)),$(round(meanB*1e3, sigdigits=4))\n")

    # close files
    map([fileX, fileY, fileZ, fileXPoints, fileYPoints, fileVariationRate]) do file
        close(file)
    end
end


# Main

minB = 1.0
maxB = 0.0
meanB = 0.0
dirName = "I=$(round(I, sigdigits=2))_N=$(round(Int, N))_z=$(round(z/h, sigdigits=2))h_y=$(round(sampleYHalfRange/h, sigdigits=2))h_x=$(round(sampleXHalfRange/l, sigdigits=2))l"

using Base.Threads
for (index, d) in enumerate(ds)
    result, minBOfZElement, maxBOfZElement, meanBOfZElement = calculateBin(planeZValue=z, distance=d)

    global minB, maxB, meanB
    minB = index == 1 ? minBOfZElement : min(minB, minBOfZElement)
    maxB = index == 1 ? maxBOfZElement : max(maxB, maxBOfZElement)
    meanB = index == 1 ? meanBOfZElement : mean((meanB, meanBOfZElement))
    variationRate = (maxB-minB)/meanB*100

    println("Current Loop Area at 2h = $(2h*100)cm, 2d = $(round(2d*100, sigdigits=3))cm, 2l = $(2l*100)cm; ")
    println("Conducting Area at x = $(round(-sampleXHalfRange/l, sigdigits=2))l~$(round(sampleXHalfRange/l, sigdigits=2))l, y = $(round(sampleYHalfRange/h, sigdigits=2))h~$(round(sampleYHalfRange/h, sigdigits=2))h; samplePoint=$(sampleIntervals):" )
    println("min B of z elment under d = $(round(d/h, sigdigits=4))h is: $(minB*1e3) [mT]")
    println("max B of z elment under d = $(round(d/h, sigdigits=4))h is: $(maxB*1e3) [mT]")
    println("Magnetic Field Variance Rate: $variationRate%\n")

    store(result; index=index, planeZValue=z, distance=d, maxB=maxB, minB=minB, meanB=meanB, variationRate=variationRate)
end