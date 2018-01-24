#include "Pattern.h"



//template <class T>
//Pattern::Pattern(T*** theBigImage, Coord startCoord, int width, int height, int depth)
//{
//	this->theBigImage = theBigImage;
//	this->startCoord = startCoord;
//	this->myWidth = width;
//	this->myHeight = height;
//	this->myDepth = depth;
//}
//
//
//template <class T>
//Pattern<T>::~Pattern()
//{
//	//  DON'T DESTROY THE_BIG_MATRIX!!!!!!!!!!!!!!!
//}
//
//
//template <class T>
//double Pattern<T>::DistanceTo(Pattern other)
//{
//	double ans = 0;										// good practice, always have a good ANS at hand.
//
//	for (int d = 0; d < depth; d++)						// now we iterate over the values of the Pattern
//	{
//		for (int h = 0; h < myHeight; h++)				// iterating...
//		{
//			for (int w = 0; w < myWidth; w++)				// still iterating...
//			{
//				// now we get the current value for THIS instance of Pattern
//				T y = theBigImage[startCoord.x + w][startCoord.y + h][startCoord.z + d];
//
//				// and now we get the current value for the instance we just received from OUTER space! (was it OTHER? yes!)
//				T x = other.theBigImage[other.startCoord.x + w][other.startCoord.y + h][other.startCoord.z + d];
//
//				// and finally, we add to ANS the square of the difference between these two misterious values!
//				ans += pow(x - y, 2);
//			}
//		}
//	}
//
//	ans = sqrt(ans);									// finally, get the square root of all evil... I meant, the ANS...
//	return ans;											// and reurn it, and go get something to eat!
//}
//
//template <class T>
//T Patern::GetValue(int width, int height, int depth)
//{
//	return theBigImage[myWidth + width][myHeight + height][myDepth + depth];
//}
