#ifndef Alignment_H
#define Alignment_H

#include <string>


class MbRandom;
class Settings;
class Tree;
class Alignment {

	public:
                                        Alignment(Settings* sp, Tree* tp, MbRandom* rp);
                                       ~Alignment(void);
                                 void   print(void);
                                 void   printToFile(std::string fn);

	private:
                                 char   convertToNuc(int x);
                                 void   setRateMatrix(double q[4][4], int whichPart);
                                 void   simulate(void);
                            MbRandom*   ranPtr;
                            Settings*   settingsPtr;
                                Tree*   treePtr;
                                int**   matrix;
                                  int   totalLength;
};

#endif