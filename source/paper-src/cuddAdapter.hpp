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

#ifndef __myKconf__cuddAdapter__
#define __myKconf__cuddAdapter__

#include "BDDAdapter.hpp"
#include "gmpxx.h"

#include <limits>
#include <cuddObj.hh>
#include <cudd.h>
#include <dddmp.h>
#include <iomanip>
#include <list>
#include <map>
#include <set>

class cuddAdapter : public BDDAdapter {
    public :
    
    cuddAdapter();
    cuddAdapter(double cacheMultiplier);
    int                 getNumVars();
    DdNode*             getZero();
    DdNode*             getOne();
    cudd::Cudd            getCudd() { return mgr; };
    const   int         nodecount();
    DdNode*             getBDD();
    int                 getLevel(DdNode* node);
    std::string         getVarName(int x);
    void                init();
    void                newVar(std::string var, std::string type);
    int                 varAtPos(int pos);

    void    readBDD(const std::string& b, const std::string &s);
    static  int MAXVAR;
    
    cudd::Cudd    mgr;
    
private:
    
    mpz_class solutions_rec(DdNode* node);
    mpf_class getProb(DdNode *node, bool parity);
    void                        shuffle(const std::vector<std::string>& order);
   
    cudd::BDD     auxReadBDD(const std::string& base, const std::string &suffix);
    void changeUp(int pos, int len, int ind, int lenhigh, int* perm);
    void changeDown(int pos, int len, int ind, int lenlow, int* perm);
    void internalRefs();
    void readVarsAndComps(const std::string& filename); 
    std::map<std::string, int>  internalReferenced, internalUsed;
    
    std::set<std::string>   quantify;
    int nroots;   
    std::vector<int>        int2extLevel, levelnodes;
    std::set<std::string>   setImplies;
    int                     countVar = 0;
    cudd::BDD               theBDD;
    void                    printVars();
    cudd::BDD               getVar(std::string var);
    std::map<std::string, std::pair<cudd::BDD, cudd::BDD> >     vars;
    std::map<std::string, std::pair<int, int> >     positions;
    std::vector<std::string>            inVars;
    std::vector<bool>                   startComponent;
    int                                 reorderMin, reorderMax, numComponents;
    std::string minVar, maxVar;
    bool                                emptyTree;
    bool                                withComponents;
    std::pair<int, int>                 findSmallestBlock(int pos, int pos2);
    int                                 numVars;
    std::set<std::pair<int, int> >      currentBlocks, presentBlocks;
    
    std::vector<mpf_class>              commonalities;
    std::vector<mpz_class>              productsPerVar;
    std::map<DdNode*, mpz_class>        tProducts, eProducts;
    std::map<DdNode*, bool>             mark;
    std::set<std::string>               artVars;
    std::map<std::string, cudd::BDD>    storage;
    
};

#endif /* defined(__myKconf__cuddAdapter__) */
