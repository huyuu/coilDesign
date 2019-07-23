import UIKit
import PlaygroundSupport
/*:
 ## Interpreted System Constants
 
*/

let pi = Double.pi

/*:
 ## Constants
 
*/

let I = 1.0
let myu0 = 4*pi*pow(10, -7)
/// 10cm
let height = 0.1
let halfLength = 0.15
let length = 2*halfLength

/*:
 ## Variables
 
*/

let points = 1000
let xSpacing = length/Double(points+1)
let ySpacing = height/Double(points+1)
let xArray = (1...points).map { (index) -> Double in
    return -halfLength + xSpacing * Double(index)
}
let yArray = (1...points).map { (index) -> Double in
    return 0 + ySpacing * Double(index)
}

/*:
 ## Functions
 
*/

func BAtPoint(x: Double, y: Double) -> Double {
    let down = 1/y * ( atan((x-halfLength)/y) - atan((x+halfLength)/y) )
    let up = 1/(y-height) * ( atan((x-halfLength)/(y-height)) - atan((x+halfLength)/(y-height)) )
    return myu0/(4*pi) * I * (-1) * (down+up) * xSpacing * ySpacing
}


/*:
 ## Main
 
*/
let start = Date()
//
//let group = DispatchGroup()
//let runQueue = DispatchQueue(label: "calculation", qos: .default, attributes: .concurrent)
//
//var result: [[Double]] = {
//    let zeroLine = [Double](repeating: 0, count: points+1)
//    let newValue = [[Double]](repeating: zeroLine, count: points+1)
//    return newValue
//}()
//result.count
//result[0].count
//
//for (xIndex, xValue) in xArray.enumerated() {
//    group.enter()
//    runQueue.async {
//        var lineResult = [Double](repeating: 0, count: points+1)
//        for (yIndex, yValue) in yArray.enumerated() {
//            result[xIndex][yIndex] = BAtPoint(x: xValue, y: yValue)
//        }
//    }
//    group.leave()
//}
//
//group.notify(queue: .global()) {
////    print("ans = \(result)")
//    print("---------------------\n")
//    let all: Double = {
//        var newValue = 0.0
//        for row in result {
//            for item in row {
//                newValue += item
//            }
//        }
//        return newValue
//    }()
//    print("all = \(all)")
//    print("time comsumption = \(Date().timeIntervalSince(start))")
//}


var result = 0.0
for (xIndex, xValue) in xArray.enumerated() {
    for (yIndex, yValue) in yArray.enumerated() {
        result += BAtPoint(x: xValue, y: yValue)
    }
}

print("result = \(result)")
print("time comsumption = \(Date().timeIntervalSince(start))")





























