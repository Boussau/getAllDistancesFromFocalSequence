#include <Bpp/Numeric/Random/RandomTools.h>
#include <Bpp/Numeric/AutoParameter.h>

#include <Bpp/App/BppApplication.h>
#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Seq/App/SequenceApplicationTools.h>

#include <Bpp/Phyl/Io.all>
#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>


using namespace std;
using namespace bpp;



// Compilation:
//g++ -o getAllDistancesFromFocalSequence getAllDistancesFromFocalSequence.cpp -std=c++0x -I/usr/local/include  -L. -L/usr/local/lib  -lbpp-core -lbpp-seq -lbpp-phyl
// Static:
//g++ -s -o getAllDistancesFromFocalSequence getAllDistancesFromFocalSequence.cpp -std=c++0x -I/usr/local/include  -L. -L/usr/local/lib -lbpp-phyl -lbpp-seq -lbpp-core --static






int main(int args, char ** argv)
{

	BppApplication getAllDistancesFromFocalSequence(args, argv, "getAllDistancesFromFocalSequence");
	getAllDistancesFromFocalSequence.startTimer();
	
	string in = ApplicationTools::getStringParameter("input.tree.file",getAllDistancesFromFocalSequence.getParams(),"none");
	
	//Read tree:
	Newick newick(true);
	TreeTemplate < Node > *tree = dynamic_cast < TreeTemplate < Node > * > (newick.read(in));

	 
	//We have the tree:
	//- need to compute pairwise distances on the tree

	string focalSeq = ApplicationTools::getStringParameter("focal.sequence",getAllDistancesFromFocalSequence.getParams(),"none");
	
	//Computing pairwise distances on the tree
	std::vector<Node*> leaves = tree->getLeaves();
	std::vector<string> leafNames;
	bool found = false;
	Node* focalNode = 0;
	size_t focalId = 0;
	for (size_t i = 0; i < leaves.size() ; ++i) {
		leafNames.push_back( leaves[i]->getName() );
		if (leafNames[i] == focalSeq) {
		  found = true;
		  focalNode =  leaves[i];
		  focalId = i;
		}
	}
	if (found == false) {
	  std::cerr << "\n\t\tError: focal sequence ("<< focalSeq << ") is not present in the tree."<<std::endl; 
	  return(-1);
	}
	std::vector<string> names;
	std::vector<double> distances;
	std::vector<size_t> numNodes;
	for (size_t j = 0;  j < leaves.size() ; ++j) {
	  if (j != focalId) {
	  //std::cout << focalSeq << " and " << leaves[j]->getName()<<std::endl;
	    //  st << focalSeq << "\t" << leaves[j]->getName() << "\t" << TreeTemplateTools::getDistanceBetweenAnyTwoNodes ( *focalNode, *(leaves[j]) ) << std::endl;
	    names.push_back(leaves[j]->getName());
	    distances.push_back(TreeTemplateTools::getDistanceBetweenAnyTwoNodes ( *focalNode, *(leaves[j]) ));
	    numNodes.push_back(TreeTemplateTools::getPathBetweenAnyTwoNodes ( *focalNode, *(leaves[j]) ).size());
	  }
	}

	std::vector<size_t> newOrder = VectorTools::order(distances);

	string out = ApplicationTools::getStringParameter("output.file",getAllDistancesFromFocalSequence.getParams(),"none");

	filebuf fb_; //file output
	fb_.open ( out.c_str(), ios::out  );
	std::ostream* outTxt_ = new ostream ( &fb_ );
	
	///	*outTxt_ << st.str() ;

	*outTxt_ <<  "focalSequence" << "\t" << "otherSequence" <<"\tdistance\tnumNodes\n";
	for (size_t i = 0 ; i < names.size() ; ++i) {
	  *outTxt_ << focalSeq << "\t" << names[newOrder[i]] << "\t" << distances[newOrder[i]] << "\t"<< numNodes[newOrder[i]]<<std::endl;
	}

	outTxt_->flush();
		
	std::cout << "getAllDistancesFromFocalSequence's done. Bye." << std::endl;
	ApplicationTools::displayTime("Total execution time:");
	
	
	return 1;
	

}
