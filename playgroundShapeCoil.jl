# Constants

const myu0 = 4pi*1e-7
const I = 1.0  # 1[A]
const N = 1
const h = 0.05  # 5cm
const l = 0.10  # 10cm
const d = 0.05  # 5cm
const R = h


# Variables

const sourceIntervals = 100  # per part
const sampleIntervals = 300


const sampleXHalfRange = l+R
const sampleXSpacing = 2*sampleXHalfRange/sampleIntervals
const sampleXIntervals = [ -sampleXHalfRange + sampleXSpacing*i for i=1:sampleIntervals-1 ]

const sampleYHalfRange = h
const sampleYSpacing = 2*sampleYHalfRange/sampleIntervals
const sampleYIntervals = [ -sampleYHalfRange + sampleYSpacing*i for i=1:sampleIntervals-1 ]

const sampleZHalfRange = d
const sampleZSpacing = 2*sampleZHalfRange/sampleIntervals
const sampleZIntervals = [ -sampleZHalfRange + sampleZSpacing*i for i=1:sampleIntervals-1 ]


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


function BAtPoint(point::Point)::Array{Float64, 1}
    result::Array{Float64, 1} = numericalIntegrateOf(c1; lowerLimit=-l, upperLimit=l, n=sourceIntervals, point=point)
    result += numericalIntegrateOf(c2; lowerLimit=-l, upperLimit=l, n=sourceIntervals, point=point)
    result += numericalIntegrateOf(c3; lowerLimit=-pi/2, upperLimit=pi/2, n=sourceIntervals, point=point)
    result += numericalIntegrateOf(c4; lowerLimit=pi/2, upperLimit=3pi/2, n=sourceIntervals, point=point)
    return result
end


# Main

# MARK: -  Calculation
const samplePoints = sampleIntervals-1
result = map(zeros((samplePoints, samplePoints, samplePoints))) do x
    Point(x, x, x)
end

@time for (xIndex, xValue) in enumerate(sampleXIntervals), (yIndex, yValue) in enumerate(sampleYIntervals), (zIndex, zValue) in enumerate(sampleZIntervals)
    global result
    resultInArray = Float64(myu0/(4pi)*(N*I)) * BAtPoint(Point(xValue, yValue, zValue))  # = BVector at P(x, y, z)
    result[xIndex, yIndex, zIndex] = Point(resultInArray)
end

println("size of result is $(size(result))")


# MARK: - Storage
# fileX = open("xElementsOfB.csv", "w")
# fileY = open("yElementsOfB.csv", "w")
# fileZ = open("zElementsOfB.csv", "w")

fileX, fileY, fileZ = map(["x", "y", "z"]) do var
    open("$(var)ElementsOfB.csv", "w")
end

for (xIndex, xValue) in enumerate(sampleXIntervals), (yIndex, yValue) in enumerate(sampleYIntervals), (zIndex, zValue) in enumerate(sampleZIntervals)
    point = result[xIndex, yIndex, zIndex]
    write(fileX, "$(round(xValue, digits=4)),$(round(yValue, digits=4)),$(round(zValue, digits=4)),$(round(point.x, sigdigits=6))\n")
    write(fileY, "$(round(xValue, digits=4)),$(round(yValue, digits=4)),$(round(zValue, digits=4)),$(round(point.y, sigdigits=6))\n")
    write(fileZ, "$(round(xValue, digits=4)),$(round(yValue, digits=4)),$(round(zValue, digits=4)),$(round(point.z, sigdigits=6))\n")
end

map([fileX, fileY, fileZ]) do file
    close(file)
end


# MARK: - Plots
# using Plots
# pyplot()
