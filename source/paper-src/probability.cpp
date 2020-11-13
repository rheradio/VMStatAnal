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

#include <stdio.h>
#include "cuddAlgo.hpp"
#include "cuddAdapter.hpp"

int main(int argc, char** argv) {
    if (argc < 2 || argc > 5) {
        std::ostringstream ost;
        ost << " Usage: " << argv[0] << " [-fast] [ -t numthreads] <bdd file>" << std::endl;
        throw std::invalid_argument(ost.str());
    }
    cuddAdapter *adapter = new cuddAdapter(1);
    int i = 1;
    bool fast = false;
    int  nt = 0;
    while (i < argc-1) {
        if (std::string(argv[i]) == "-fast") {
            fast = true;
            i++;
            continue;
        }
        if (std::string(argv[i]) == "-t") {
            nt = atoi(argv[i+1]);
            i = i+2;
            continue;
        }
    }
    
    adapter->readBDD(argv[i], "");
    std::vector<mpf_class> com;
    
    if (nt == 0) {
        com = computeCommonality(adapter);
    }
    else {
        com = computeCommonalityMP(nt, adapter, fast);
    }
    
    int x = 0;
    signed int max = 0;
    for(int l = 0; l < adapter->getNumVars(); l++)
        if (adapter->getVarName(l).length() > max)
            max = adapter->getVarName(l).length();
    for(mpf_class& f : com)
        std::cout << std::setw(max+1) << adapter->getVarName(x++) << " " << f << std::endl;
}
    
