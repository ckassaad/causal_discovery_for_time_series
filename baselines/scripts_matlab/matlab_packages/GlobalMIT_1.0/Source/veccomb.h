#include <vector>
#include "combination.h"


template<class T>
bool FindVecComb(
	std::vector< std::vector< T > >& VecN,
	std::vector< std::vector< T > >& VecR
	)
{
	if(VecR.size() == 0)
		return false;

	size_t nDepth = VecR.size()-1;
	bool NoMore = true;
	do
	{
		NoMore = stdcomb::next_combination(	
			VecN.at(nDepth).begin(), 
			VecN.at(nDepth).end(), 
			VecR.at(nDepth).begin(), 
			VecR.at(nDepth).end() );

		if(NoMore == false)
		{
			// Initailize as before
			for( size_t i=0; i<VecR.at(nDepth).size(); ++i)
			{
				VecR.at(nDepth).at(i) = VecN.at(nDepth).at(i);
			}
		}

		if( nDepth == 0 && NoMore == false)
			return false;

		--nDepth;
	}
	while(NoMore == false);

	return true;
}
