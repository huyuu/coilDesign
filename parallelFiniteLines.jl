# Constants

const myu0 = 4pi*1e-7
const I = 1.0  # 1[A]
const height = 0.05  # 5cm
const halfLength = 0.1  # 10cm
const length = 2*halfLength


# Variables

const points = 10000
const xSpacing = length/(points+1)
const ySpacing = height/(points+1)
const xArray = [ Float64(-halfLength + xSpacing*index) for index=1:points ]
const yArray = [ Float64(ySpacing*index) for index=1:points ]


# Functions

function BAtPoint(x::Float64, y::Float64)::Float64
    down = 1/y * ( atan((x-halfLength)/y) - atan((x+halfLength)/y) )
    up = 1/(y-height) * ( atan((x-halfLength)/(y-height)) - atan((x+halfLength)/(y-height)) )
    myu0/(4*pi) * I * (-1) * (down+up) * xSpacing * ySpacing
end


# Main

@time result = [ BAtPoint(x, y) for x in xArray, y in yArray ]
println("size of result is $(size(result))")
println("sum of B is $(sum(result))")
