#ifndef Node_H
#define Node_H



class Node {

	public:
                                        Node(void);
                               double   getProportion(void) { return proportion; }
                                Node*   getLft(void) { return lft; }
                                Node*   getRht(void) { return rht; }
                                Node*   getAnc(void) { return anc; }
                                  int   getIndex(void) { return index; }
                                  int   getNuc(void) { return nuc; }
                                 void   setProportion(double x) { proportion = x; }
                                 void   setLft(Node* p) { lft = p; }
                                 void   setRht(Node* p) { rht = p; }
                                 void   setAnc(Node* p) { anc = p; }
                                 void   setIndex(int x) { index = x; }
                                 void   setNuc(int x) { nuc  = x; }

	private:
                               double   proportion;
                                Node*   lft;
                                Node*   rht;
                                Node*   anc;
                                  int   index;
                                  int   nuc;
};

#endif