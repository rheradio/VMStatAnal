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

#include "cuddAlgo.hpp"


std::map<DdNode*, mpf_class>                                        probabilities;
std::list<DdNode*>                                                  buffer;
std::vector<std::list<std::pair<short int, short int>>::iterator>   implierUp, implierDown, excluder;
std::list<std::pair<short int, short int>>                          implicationsUp, implicationsDown, exclusions;

void getProducts(int plevel, int tlevel, int elevel, const mpz_class& tr, const mpz_class& er, mpz_class& thenPart, mpz_class& elsePart) {
    // thenPart = 2 ^ (tlevel-plevel-1) * tr
    mpz_ui_pow_ui(thenPart.get_mpz_t(), 2, tlevel - plevel -1);
    thenPart *= tr;
    
    // elsePart = 2 ^(elevel - plevel -1) * er
    mpz_ui_pow_ui(elsePart.get_mpz_t(), 2, elevel - plevel -1);
    elsePart *= er;
}

mpz_class _count(int plevel, int tlevel, int elevel, const mpz_class& tr, const mpz_class& er) {
    mpz_class thenPart, elsePart;
    getProducts(plevel, tlevel, elevel, tr, er, thenPart, elsePart);
    
    return (thenPart + elsePart);
}

mpz_class count(cuddAdapter* a) {
    return traverse(a, mpz_class(0), mpz_class(1), &_count);
}

int _listNodes(int plevel, int tlevel, int elevel, const int& tr, const int& er, DdNode* node) {
    buffer.push_back(node);
    return 0;
}
mpz_class countMP(int threads, cuddAdapter* a, bool fast) {
    return traverseMP(threads, a, mpz_class(0), mpz_class(1), &_count, fast);
}
mpz_class _probabilities(int plevel, int tlevel, int elevel, const mpz_class& tr, const mpz_class& er, DdNode* node) {
    mpz_class thenPart, elsePart;
    getProducts(plevel, tlevel, elevel, tr, er, thenPart, elsePart);
    mpf_class fthen = thenPart;
    mpf_class felse = elsePart;
    probabilities[node] = fthen / (fthen + felse);
    return (thenPart + elsePart);
}

void compProbabilitiesMP(int threads, cuddAdapter* a, bool fast) {
    traverseMP(threads, a, mpz_class(0), mpz_class(1), &_probabilities, fast);
}

void compProbabilities(cuddAdapter* a) {
    traverse(a, mpz_class(0), mpz_class(1), &_probabilities);
}


int distcounter = 1;
int totalNodes  = 0;
int divider     = 0;
std::vector<std::vector<mpz_class>> combinations;
std::vector<mpz_class> powersof2;
std::set<int> levelJump;

mpz_class& getPowerof2(int n) {
  return powersof2[n];
}

void makePowersof2(int n ) {
  powersof2.resize(n+1);
  levelJump.insert(0);
  for(int i : levelJump)
     mpz_ui_pow_ui(powersof2[i].get_mpz_t(), 2, i);
}
 
mpz_class& getCombinations(int a, int b) {
    if (b < a/2+1) {
        return combinations.at(a).at(b);
    }
    else {
        return combinations.at(a).at(a-b);
    }
}
void makeCombinations(int n) {
    std::vector<mpz_class> lastRow, current;
    mpz_class lastNum;
    // 0 over 0
    lastRow.push_back(1);
    combinations.push_back(lastRow);
    lastRow.clear();
    // 1 over 0
    lastRow.push_back(1);
    // 1 over 1
    lastRow.push_back(1);
    combinations.push_back(lastRow);
    
    for(int x = 2; x <= n; x++) {
        if (x % 1000 == 0) std::cerr << "combinations: " << x << std::endl;
        current.clear();
        lastNum = 0;
        for(int k = 0; k < x/2+1; k++) {
            mpz_class trav = getCombinations(x-1, k);
            current.push_back(lastNum+trav);
            lastNum = trav;
        }
        combinations.push_back(current);
        lastRow = current;
    }
}
//void makeCombinations(int n) {
//    mpz_class factor   = 1;
//    int multiply       = n;
//    int divide         = 1;
//    int half           = n/2;
//    for(int x = 0; x <= half; x++){
//        combMap[n].push_back(factor);
//        factor = (factor*multiply--)/divide++;
//    }
//    if (n % 2 == 0) half--;
//    for(int x = half;  x >= 0; x--) {
//        combMap[n].push_back(combMap[n][x]);
//    }
//}

// We return the biggest jump.
int  _compLevelJump(int plevel, int tlevel, int elevel,
                                     const int& tr,
                                     const int& er) {
    int max = -1;
    if (tr > max) max = tr;
    if (er > max) max = er;

    if (tr != -1) {
	levelJump.insert(tlevel - plevel - 1);
        if (tlevel - plevel - 1 > max) {
	    max = tlevel - plevel - 1;
        }
    }
    if (er != -1) {
        levelJump.insert(elevel - plevel - 1);
        if (elevel - plevel - 1 > max) {
	    max = elevel - plevel - 1;
       }
    }
    if (max == -1) max = 0;
    return max;
}
    
std::vector<mpz_class>  _distCombine(int plevel, int tlevel, int elevel,
                                     const std::vector<mpz_class>& tr,
                                     const std::vector<mpz_class>& er) {
    if (distcounter % divider == 0) 
       if (1+100*(distcounter)/totalNodes <= 100)
           std::cerr << 1+100*(distcounter)/totalNodes << "%" << std::endl;
    distcounter++;
    std::vector<mpz_class> thenBit(tlevel - plevel + tr.size());
    std::vector<mpz_class> elseBit(elevel - plevel + er.size());
    // Size of the distribution is any of the lowdif or highdif
    // it's the same number
    std::vector<mpz_class> infoDist(elseBit.size());
    // Now we compute the distribution of the "else" (low) path.
    
    // If the else node is zero, then we don't do anything
    if (er.size() != 1 || er[0] != 0)
        for(int a = 0; a < elseBit.size(); a++) {
            // We fill theBit one component at a time
            for(int b = 0;  b < elevel - plevel ; b++) {
                // How many different paths (along omitted nodes) are there from
                // the node to the "else" descendant passing through b omitted
                //  "then" arcs? Answer ((elevel - plevel -1) over b)
                // factor is (elevel - plevel -1) over b
                if (a-b >= 0 && a-b < er.size()) {
                    elseBit[a] += getCombinations(elevel-plevel-1, b) * er[a-b];
                }
            }
        }
    
    // Now we compute the distribution of the "then" (high) path. (almost symmetrical)
    // a starts at 1 because thenBit[0] has to be 0.
    // If the then node is zero, then we don't do anything
    if (tr.size() != 1 || tr[0] != 0)
        for(int a = 1; a < thenBit.size(); a++) {
            // We fill highdif one component at a time
            for(int b = 0; b < tlevel - plevel; b++) {
                // How many different paths (along omitted nodes) are there from
                // the node to the "then" descendant passing through b omitted
                //  "then" arcs? Answer (tlevel - plevel - 1) over b
                // factor is (tlevel - plevel - 1) over b
                if (a-b >= 0 && a-b-1 < tr.size()) {
                    thenBit[a] += getCombinations(tlevel-plevel-1, b) * tr[a-b-1];
                }
            }
        }
    // Now we combine both;
    for(int l = 0; l < infoDist.size(); l++) {
        // We add the branches
        infoDist[l] = thenBit[l] + elseBit[l];
    }
    return infoDist;
}


std::vector<mpz_class>  computeDistributionMP(int threads, cuddAdapter* a, bool fast) {
    std::vector<mpz_class> dZero, dOne;
    dZero.resize(1);
    dOne.resize(1);
    dZero[0] = 0;
    dOne[0]  = 1;
    // Jumps to node zero are ignored. But jumps to one are counted.
    int maxJump = traverse(a, -1, 0, _compLevelJump);
    totalNodes = a->nodecount();
    divider    = totalNodes/100;
    makeCombinations(maxJump);
    return traverseMP(threads, a, dZero, dOne, _distCombine, fast);
}
std::vector<mpz_class> computeDistribution(cuddAdapter* a) {
    std::vector<mpz_class> dZero, dOne;
    dZero.resize(1);
    dOne.resize(1);
    dZero[0] = 0;
    dOne[0]  = 1;
    int maxJump = traverse(a, -1, 0, _compLevelJump);
    totalNodes = a->nodecount();
    divider    = totalNodes/100;
    makeCombinations(maxJump);
    return traverse(a, dZero, dOne, _distCombine);
}

class ucr {
public:
    // How many solutions from this node
    mpz_class               counter;
    // For each variable from this level to the bottom, how many
    // solutions with the variable set to true
    std::vector<mpz_class>  products;
    
};
int counter  = 1;
ucr _upProducts(int plevel, int tlevel, int elevel, const ucr& tr, const ucr& er) {
    ucr res, resT, resE;
    mpz_class thenPart, elsePart;
    if (counter % divider == 0)
        if (1+100*(counter)/totalNodes <= 100)
            std::cerr << 1+100*(counter)/totalNodes << "%" << std::endl;
    counter++;
    if (tr.counter != 0) {
        mpz_class& tempThen = getPowerof2(tlevel - plevel - 1);
        thenPart = tr.counter*tempThen;
        mpz_class halfThen = thenPart/2;
	resT.counter = thenPart;
        resT.products.push_back(thenPart);
        // For reduced nodes in to the then child
        for(int i = plevel+1; i < tlevel; i++)
            resT.products.push_back(halfThen);

        // temp = 2 ^ (tlevel - plevel -1)
        for(const mpz_class& t : tr.products) {
            resT.products.push_back(t*tempThen);
        }
    }
    
    if (er.counter != 0) {
        mpz_class& tempElse = getPowerof2(elevel - plevel - 1);
        elsePart = er.counter*tempElse;
        mpz_class halfElse = elsePart/2;
	resE.counter = elsePart;
        resE.products.push_back(0);
        for(int i = plevel+1; i < elevel; i++)
            resE.products.push_back(halfElse);

        for(const mpz_class& e : er.products) {
            resE.products.push_back(e*tempElse);
        }
    }
    if (thenPart != 0 && elsePart != 0) {
        res.counter = thenPart + elsePart;
	for(int i = 0; i < resT.products.size(); i++) 
           res.products.push_back(resT.products[i] + resE.products[i]);

       return res;
   }
   if (thenPart == 0) 
	return resE;
   
  return resT;   
}



std::vector<mpf_class> computeCommonality(cuddAdapter *a) {
    std::vector<mpf_class> res;
    ucr rZero, rOne, recRes;
    rZero.counter = 0;
    rOne.counter  = 1;
    rZero.products.push_back(0);
    rOne.products.push_back(1);
    // Jumps to node zero are ignored. But jumps to one are counted.
    int maxJump = traverse(a, -1, 0, _compLevelJump);
    makePowersof2(maxJump);
    totalNodes = a->nodecount();
    divider    = totalNodes/100;
    recRes = traverse(a, rZero, rOne, _upProducts);
    // Last component would be commonality of 1 so we do not include it
    for(int i = 0; i < recRes.products.size()-1; i++) {
        res.push_back((mpf_class) recRes.products[i] /recRes.counter);
}
    return res;
}

std::vector<mpf_class> computeCommonalityMP(int num, cuddAdapter *a, bool fast) {
    std::vector<mpf_class> res;
    ucr rZero, rOne, recRes;
    rZero.counter = 0;
    rOne.counter  = 1;
    rZero.products.push_back(0);
    rOne.products.push_back(1);
    // Jumps to node zero are ignored. But jumps to one are counted.
    int maxJump = traverse(a, -1, 0, _compLevelJump);
    makePowersof2(maxJump);
    totalNodes = a->nodecount();
    divider    = totalNodes/100;
    recRes = traverseMP(num, a, rZero, rOne, _upProducts, fast);
    // Last component would be commonality of 1 so we do not include it
    for(int i = 0; i < recRes.products.size()-1; i++)
        res.push_back((mpf_class) recRes.products[i] /recRes.counter);
    
    return res;
}

std::vector<int> levelNodes;

int _combLevelNodes(int plevel, int tlevel, int elevel,
                                 const int& tr,
                                 const int& er) {
    levelNodes[plevel]++;
    return 0;
}

std::vector<int> getLevelNodes(cuddAdapter* a) {
    levelNodes.resize(a->getNumVars());
    for(int i = 0; i < a->getNumVars(); i++)
        levelNodes[i] = 0;
    traverse(a, 0, 0, _combLevelNodes, false);
    return levelNodes;
}

