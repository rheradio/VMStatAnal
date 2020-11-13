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

#include "cuddAdapter.hpp"
#include "cuddInt.h"
#include "cuddAlgo.hpp"

int cuddAdapter::MAXVAR = 100000;

DdNode*          cuddAdapter::getBDD()  { return theBDD.getNode();         }
DdNode*          cuddAdapter::getZero() { return mgr.bddZero().getNode();  }
DdNode*          cuddAdapter::getOne()  { return mgr.bddOne().getNode();   }
std::vector<int> cuddAdapter::pos2var() {
    std::vector<int> res;
    for(int x = 0; x < numVars; x++)
        res.push_back(mgr.ReadInvPerm(x));
    return res;
}
int    cuddAdapter::varAtPos(int pos) {
    return mgr.ReadInvPerm(pos);
}
std::string        cuddAdapter::getVarName(int x) { return inVars[x]; }
int                 cuddAdapter::getNumVars()        { return numVars;                   }
void cuddAdapter::init() {
    numVars = 0;
    theBDD = mgr.bddOne();
    minVar = maxVar = "";
}
int cuddAdapter::getLevel(DdNode* node) {
    if (Cudd_IsConstant(Cudd_Regular(node))) {
        return (mgr.ReadSize());
    }
    else {
        return (mgr.ReadPerm(Cudd_Regular(node)->index));
    }
}

void   cuddAdapter::newVar(std::string var, std::string type) {
    //pcomp->newVariable();
    countVar++;
    vars[var] = std::pair<cudd::BDD, cudd::BDD>(mgr.bddVar(), mgr.bddZero());
    positions[var] = std::pair<int, int>(numVars, 0);
    inVars.push_back(var);
    mgr.pushVariableName(var);
    numVars++;
}
cuddAdapter::cuddAdapter(double cacheMultiplier) :
        mgr(0,0,CUDD_UNIQUE_SLOTS, cacheMultiplier*CUDD_CACHE_SLOTS, 0) {
            withComponents = false;
            init();
}

cuddAdapter::cuddAdapter() :
mgr(0,0,CUDD_UNIQUE_SLOTS, CUDD_CACHE_SLOTS, 0) {
    withComponents = false;
}

const int cuddAdapter::nodecount() {
    return Cudd_DagSize(theBDD.getNode());
}
void cuddAdapter::readVarsAndComps(const std::string& filename) {
    // The variables must be created in the manager before we read the file
    // So we scan the header for the info we need and create them in case
    // they don't exist yet.
    
    if (mgr.ReadSize() == 0) {
        int nvars = 0;
        // First we read the header
        std::ifstream dump(filename);
        if (!dump.good()) {
            throw std::invalid_argument("Couldn't open file "+filename);
        }
        
        std::string word;
        std::vector<std::string> shuf;
        for(;;) {
            if (dump >> word) {
                if (word == ".nroots") {
                    dump >> nroots;
                }
                else
                if (word == ".nvars") {
                    dump >> nvars;
                }
                // Now it's not support variable names, it's ALL variable names
                else if (word == ".varnames")
                    for(int x = 0; x < nvars; x++) {
                        dump >> word;
                        newVar(word, "boolean");
                    }
                else if (word == ".orderedvarnames") {
                          for(int x = 0; x < nvars; x++) {
                               dump >> word;
                              shuf.push_back(word);
                          }
                    shuffle(shuf);
                }
                else if (word == ".nodes")
                    break;
            }
        }
        dump.close();
    }
}

cudd::BDD   cuddAdapter::auxReadBDD(const std::string& base, const std::string &suffix)
{
    std::string filename;
    if (suffix == "")
        filename = base + ".dddmp";
    else
        filename = base + "-" + suffix + ".dddmp";
    
    readVarsAndComps(filename);
    Dddmp_VarMatchType  varMatchMode = DDDMP_VAR_MATCHPERMIDS;
    FILE *fp = fopen(filename.c_str(), "r");
    int mode = DDDMP_MODE_TEXT;
    return cudd::BDD(mgr,
                     Dddmp_cuddBddLoad(mgr.getManager(),
                                       varMatchMode,
                                       NULL,
                                       NULL,
                                       NULL,
                                       mode,
                                       (char*)filename.c_str(),
                                       fp)
                     );
}

void    cuddAdapter::readBDD(const std::string& base, const std::string &suffix) {
    theBDD = auxReadBDD(base, suffix);
}

void    cuddAdapter::shuffle(const std::vector<std::string>& extOrder) {
    // If there are no variables...
    if (extOrder.empty())
        return;
    if (extOrder.size() != numVars) {
        std::ostringstream ost;
        ost << "Shuffle size " << extOrder.size() << " != " << numVars << " numVars" << std::endl;
        if (numVars == 0) ost << "Maybe you are reading a dddmp file in an old format (without variable names)";
        for(const std::string& s : extOrder) std::cerr << "extOrder " << s << std::endl;
        throw std::logic_error(ost.str());
    }
    int intOrder[numVars], cont = 0;
    std::set<int> check;
    for(const std::string& s : extOrder)
        if (positions[s].second == 0) {
            intOrder[cont++] = positions[s].first;
        }
        else {
            intOrder[cont++] = positions[s].first;
            intOrder[cont++] = positions[s].second;
        }
    //mgr.SetTree(Mtr_InitGroupTree(mgr.ReadInvPerm(0), numVars));
    mgr.MyShuffleHeap(intOrder);
};
