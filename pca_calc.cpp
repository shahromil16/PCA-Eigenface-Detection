#include <shark/Algorithms/Trainers/PCA.h>
#include <shark/Data/Pgm.h>
#include <fstream>
#include <boost/format.hpp>
#include <math.h> 
#include <matrix.h> 
#include <mex.h>

using namespace std;
using namespace shark;


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) - See more at: http://www.shawnlankton.com/2008/03/getting-started-with-mex-a-short-tutorial/#sthash.BK94hYqM.dpuf
{
	mxArray *images = prhs[0];
	float a, b, c;
	unsigned l = images.numberOfElements();   // number of samples
	unsigned x = imagesInfo.element(0).x; // width of images
	unsigned y = imagesInfo.element(0).y; // height of images

	cout << "Eigenvalue" << flush;
	PCA pca(images);
	
	cout << "Mean and eigenfaces" << flush;
	ofstream ofs("facesEigenvalues.csv");
	for (unsigned i = 0; i < l; i++)
		ofs << pca.eigenvalue(i) << endl;
	
	unsigned m = 200;
	LinearModel<> enc;
	pca.encoder(enc, m);
	Data<RealVector> encodedImages = enc(images);
	
	LnearModel<> dec;
	pca.decoder(dec, m);
	boost::format fmterRec("Reconstructedfaces.pgm");
	exportPGM((fmterRec % sImg % m).str().c_str(), dec(encodedImages.element(sImg)), x, y);
	for (unsigned i = 0; i < l; i++)
	a << pca.eigenvalue(i);
	b << pca.eigenvectors(i); << endl;

	plhs[0] = mxCreateDoubleMatrix(a)
	plhs[1] = mxCreareDoubleMatrix(b)
	plhs[2] = mxCreateDoubleMatrix(c)

}