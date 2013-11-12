import Control.Monad
import Data.List

data F3 = Zero | One | Two deriving (Show)

data Matrix x = Matrix { a :: x, b :: x, c :: x, d :: x } deriving (Eq)

identity = Matrix One Zero Zero One
 
multiplyInField  Zero _ = Zero
multiplyInField  _ Zero = Zero
multiplyInField One Two = Two
multiplyInField Two One = Two
multiplyInField One One = One
multiplyInField Two Two = One

addInField  Zero a = a
addInField  a Zero = a
addInField One Two = Zero
addInField Two One = Zero
addInField One One = Two
addInField Two Two = One

multiplicativeInverse One = One
multiplicativeInverse Two = Two

instance Eq F3 where
	Zero == Zero = True
	Zero == One = False
	One == Zero = False
	One == Two = False
	Two == One = False
	Zero == Two = False
	Two == Zero = False
	One == One = True
	Two == Two = True

instance (Show a) => Show (Matrix a) where
	show (Matrix a b c d) = "[ " ++ show a ++ " " ++ show b ++ " ]\n[ " ++ show c ++ " " ++ show d ++ " ]\n"

instance Num F3 where
	(+) = addInField
	(*) = multiplyInField

	abs = id

	negate Zero = Zero
	negate One = Two
	negate Two = One

	signum _ = One

	fromInteger a
		| a `mod` 3 == 0 = Zero
		| a `mod` 3 == 1 = One
		| a `mod` 3 == 2 = Two

-- matrix multiplication
multiply :: (Num a) => Matrix a -> Matrix a -> Matrix a
multiply one another = Matrix (a1*a2 + b1*c2) (a1*b2 + b1*d2) (c1*a2 + d1*c2) (c1*b2 + d1*d2) where
	a1 = a one
	b1 = b one
	c1 = c one
	d1 = d one
	a2 = a another
	b2 = b another
	c2 = c another
	d2 = d another

det :: (Num a) => Matrix a -> a
det m = (a m)*(d m) - (b m)*(c m)

toMatrix :: [a] -> Matrix a
toMatrix x = Matrix (x !! 0) (x !! 1) (x !! 2) (x !! 3)

isInvertible :: Matrix F3 -> Bool
isInvertible m = (det m) /= Zero

isSpecial :: Matrix F3 -> Bool
isSpecial m = det m == One

-- Use the well-known formula for the inverse of an 2*2 matrix, using that all of the matrices in the
-- special linear group have determinant 1, and so 1/det(m) = 1.
getInverseOfSpecialMatrix :: (Num a) => Matrix a -> Matrix a
getInverseOfSpecialMatrix (Matrix a b c d) = Matrix d (-b) (-c) a

conjugate :: (Num a) => Matrix a -> Matrix a -> Matrix a
conjugate x m = (m `multiply` x) `multiply` (getInverseOfSpecialMatrix m)

-- We generate the conjugacy class for a single element by conjugating it by all elements in SL_2(F_3)
-- and then removing the duplicate elements.  We implicitly use the fact that conjugating an element by every element of
-- the group will give us its entire conjugacy class. 
conjugacyClass :: (Num a, Eq a) => [Matrix a] -> Matrix a -> [Matrix a]
conjugacyClass group m = nub $ map (conjugate m) group

listsAreIdentical :: (Eq a) => [a] -> [a] -> Bool
listsAreIdentical x y 
	| length x /= length y = False
	| otherwise = all (True==) $ map (\a -> a `elem` y) x

{- Matices and lists of matrices -}

specialLinearGroupOnF3 :: [Matrix F3]
specialLinearGroupOnF3 = filter isSpecial $ map toMatrix $ replicateM 4 [Zero, One, Two]

specialLinearInverses :: [Matrix F3]
specialLinearInverses = map getInverseOfSpecialMatrix specialLinearGroupOnF3

-- We generate the set of conjugacy classes for each element in SL_2(F_3) and then remove the duplicate conjugacy classes
-- to render a list of distinct conjugacy classes.
specialLinearConjugacyClasses :: [[ Matrix F3 ]]
specialLinearConjugacyClasses = nubBy listsAreIdentical $ map (conjugacyClass specialLinearGroupOnF3) specialLinearGroupOnF3
