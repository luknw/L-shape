import System.Environment

import Numeric.LinearAlgebra
import Graphics.EasyPlot

localMatrix :: Int -> Matrix Double
localMatrix n =
    let
        scaleFactor = (1 / elementEdgeLen n)^2
        unscaled = fromLists [ [2/3, -1/6, -1/3, -1/6]
                             , [-1/6, 2/3, -1/6, -1/3]
                             , [-1/3, -1/6, 2/3, -1/6]
                             , [-1/6, -1/3, -1/6, 2/3] ]
    in
        scale scaleFactor unscaled

element2Indices :: Int -> Int -> [((Int,Int), Double)]
element2Indices n position = zip (elementIndicesSquared n position) (toList . flatten $ localMatrix n)
  where
    leftTopIndex n p
        | p <= 2 * n^2 = p `div` (2*n + 1) + p - 1
        | otherwise    = (p - (2 * n^2 + 1)) `div` n + p + 2*n - 1
    increment n p
        | p <= 2 * n^2 = 2*n + 1
        | otherwise    = n + 1
    elementIndices n p = [ leftTopIndex n p
                         , (leftTopIndex n p) + 1
                         , (leftTopIndex n p) + (increment n p) + 1
                         , (leftTopIndex n p) + (increment n p)]
    elementIndicesSquared n p = [(i, j) | i <- elementIndices n p, j <- elementIndices n p]

considerDirichletBoundaries :: Int -> (Matrix Double, Vector Double) -> (Matrix Double, Vector Double)
considerDirichletBoundaries n (m,v) =
    let
        row = n + 1
        horizontalBoundary = take row [(2*n + 1)*n ..]
        commonPoint = last horizontalBoundary
        verticalBoundary = take n [commonPoint + row, commonPoint + 2*row ..]
        boundaryIndices = horizontalBoundary ++ verticalBoundary

        baseVector :: Int -> Int -> Vector Double
        baseVector dims i = assoc dims 0 [(i, 1)]

        replaceRow :: Int -> Matrix Double -> Vector Double -> Matrix Double
        replaceRow i m new = (m ? [0 .. i-1]) === asRow new === (m ? [i+1 .. rows m - 1])

        replaceWithBaseRow i m = replaceRow i m (baseVector (cols m) i)

        m' = foldr replaceWithBaseRow (m ||| asColumn v) boundaryIndices
    in
        (takeColumns (cols m' - 1) m', flatten $ dropColumns (cols m' - 1) m')

globalPointsCount :: Int -> Int
globalPointsCount n = (n + 1) * (3*n + 1)

globalMatrix :: Int -> Matrix Double
globalMatrix n =
    let
        elementCount = 3 * n^2
        globalPoints = concat $ map (element2Indices n) [1 .. elementCount]

        zeros = (globalPointsCount n><globalPointsCount n) $ repeat 0
    in
        accum zeros (+) globalPoints

globalIndex2Position :: Int -> Int -> (Double, Double)
globalIndex2Position n i
    | i <= lastUpperPointIndex n =
        let
            (d, m) = i `divMod` (2*n + 1)
         in
            ( (-1) + (fromIntegral m) * (elementEdgeLen n)
            , 1    - (fromIntegral d) * (elementEdgeLen n))
    | otherwise =
        let
            (d, m) = (i - (2*n + 1) * (n+1)) `divMod` (n+1)
         in
            ( 0                   + (fromIntegral m) * (elementEdgeLen n)
            , (-elementEdgeLen n) - (fromIntegral d) * (elementEdgeLen n))

g :: (Double, Double) -> Double
g (x,y) = x

lastUpperPointIndex :: Int -> Int
lastUpperPointIndex n = (2*n + 1)*(n+1) - 1

neumannBoundaryIndices :: Int -> [Int]
neumannBoundaryIndices n =
    let
        lui = lastUpperPointIndex n

        upRowCount i
            | i <= lui  = 2*n + 1
            | otherwise = n + 1

        up i = i - upRowCount i

        downRowCount i
            | i <= up lui   = upRowCount i
            | otherwise     = upRowCount (lui + 1)

        down i = i + downRowCount i

        downFrom i = down i : map down (downFrom i)
    in
        take n [0, down 0 ..]
        ++ take (2*n) [1..]
        ++ take (2*n) (downFrom (2*n))
        ++ take (n-1) [globalPointsCount n - 2, globalPointsCount n - 3 ..]

elementEdgeLen :: Int -> Double
elementEdgeLen n = recip (fromIntegral n)

edgePositions :: Int -> Int -> [(Double, Double)]
edgePositions n i
    | i == 0           = [(-1, 1 - halfEdge n), (-1 + halfEdge n, 1)]
    | i == 2*n         = [(1 - halfEdge n, 1), (1, 1 - halfEdge n)]
    | i == lastIndex n = [(1, -1 + halfEdge n), (1 - halfEdge n, -1)]
    | otherwise =
        let
            x = fst $ globalIndex2Position n i
            y = snd $ globalIndex2Position n i

            eps = 1e-6
            d a b = abs (a - b)
            (~=) a b = d a b < eps

            comparePosition a b
                | x ~= (-1) = [(-1, y - halfEdge n), (-1, y + halfEdge n)]
                | x ~= 1    = [(1, y + halfEdge n), (1, y - halfEdge n)]
                | y ~= (-1) = [(x + halfEdge n, -1), (x - halfEdge n, -1)]
                | y ~= 1    = [(x - halfEdge n, 1), (x + halfEdge n, 1)]
        in
            comparePosition x y
  where
    lastIndex n = globalPointsCount n - 1
    halfEdge n = elementEdgeLen n / 2

neumannIntegral :: Int -> Int -> Double
neumannIntegral n i =
    (0.5 * elementEdgeLen n *) . foldr1 (+) . map g $ edgePositions n i
    
rhsVector :: Int -> Vector Double
rhsVector n =
    assoc (globalPointsCount n) 0 [(i, neumannIntegral n i) | i <- neumannBoundaryIndices n]

solve :: Int -> Vector Double
solve n =
    let
        m = globalMatrix n
        v = rhsVector n
        (m',v') = considerDirichletBoundaries n (m, v)
    in
        m' <\> v'

main' :: Int -> IO ()
main' n =
    let
        xys = map (globalIndex2Position n) [0 .. (globalPointsCount n) - 1]
        zs = toList $ solve n
        xyzs = zipWith (\(x,y) z -> (x,y,z)) xys zs
    in
        (plot X11 $ xyzs) >>= \b
        -> return ()

main :: IO ()
main = do
    (n:_) <- getArgs
    main' (read n :: Int)
