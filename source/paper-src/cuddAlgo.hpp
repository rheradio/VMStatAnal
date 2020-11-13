/* This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

// Written by David Fernandez Amoros 2019

#ifndef cuddAlgo_hpp
#define cuddAlgo_hpp

#include <stdio.h>
#include <sstream>
#include "gmp.h"

#include "cuddAdapter.hpp"
#include "Traverser.hpp"

mpz_class                         count(cuddAdapter* f);
mpz_class                         countMP(int threads, cuddAdapter* f, bool fast);

void                              compProbabilities(cuddAdapter* a);
void                              compProbabilitiesMP(int threads, cuddAdapter* a, bool fast);


mpf_class                         getCommonality(int i);
std::vector<mpf_class>            computeCommonality(cuddAdapter *a);
std::vector<mpf_class>            computeCommonalityMP(int num, cuddAdapter *a, bool fast);

void                              commonality(cuddAdapter* a);

std::vector<mpz_class>            computeDistribution(cuddAdapter* a);
std::vector<mpz_class>            computeDistributionMP(int threads, cuddAdapter* a, bool fast);
std::vector<mpz_class>            crappyDistribution(cuddAdapter* a);


std::vector<int>                  getLevelNodes(cuddAdapter* a);


#endif /* cuddAlgo_hpp */
