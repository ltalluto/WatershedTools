#include <Rcpp.h>
using namespace Rcpp;

NumericVector multiply_sparse(NumericMatrix mat, NumericVector y);

//' Compute a concentration derivative for an entire river network
//' @param t the time step
//' @param y the state vector
//' @param adjacencyQ A 3-column matrix giving the rows and columns of the non-zero entries of 
//'     the adjacency matrix (columns 1 and 2) and the discharge value (column 3)
//' @param qout A vector of output discharge
//' @param qin A vector of input discharge
//' @param lateral A vector giving concentraion of lateral input
//' @param csArea Cross Sectional area of each site
//' @param dx Length of each site
// [[Rcpp::export]]
NumericVector dCdt_transport_cpp(double t, NumericVector y, NumericMatrix adjacencyQ, 
		NumericVector qout, NumericVector qin, NumericVector lateral, NumericVector csArea, 
		NumericVector dx)
{
	NumericVector inputMass = multiply_sparse(adjacencyQ, y);
	NumericVector totalInputMass = inputMass + (qout - qin) * lateral;
	NumericVector totalOutputMass = qout * y;
	NumericVector advection = (-1/csArea) * (totalOutputMass - totalInputMass)/dx;
	return advection;
}



//' Quick and dirty vector by sparse matrix multiplication
//' note that no bounds checking is done
//' @name multiply_sparse
//' @param mat A 3-column matrix giving the rows and columns of the non-zero entries of the 
//'     sparse matrix (columns 1 and 2) and the value (column 3)
//' @param y The vector by which to multiply mat
//' @keywords internal
NumericVector multiply_sparse(NumericMatrix mat, NumericVector y)
{
	NumericVector output(y.size(), 0.0);
	for(size_t i = 0; i < mat.nrow(); i++) {
		// subtracting 1 because C++ is zero indexed and these indices come from R
		size_t j = mat(i,0) - 1; 
		size_t k = mat(i,1) - 1;
		output[j] += y[k] * mat(i,2);
	}
	return output;
}

//' Connect points in a watershed
//' @param dsPixel A vector of downstream pixels, diPixel[i] is downstream from pixel i
// [[Rcpp::export]]
std::vector<int> connectCPP(NumericVector dsPixel, int upstream, int downstream) {
	std::vector<int> connectedPts;
	int current = upstream;
	connectedPts.push_back(current);
	while(current != downstream && current > 1) {
		current = dsPixel(current - 1); // -1 to correct for C++ indexing
		connectedPts.push_back(current);
	}
	if(current == 1 && downstream != 1)
		connectedPts.clear();
	return connectedPts;
}


//' Construct a distance matrix between one set of points and another
//'
//' @param x  Set of pixels to which to compute distance
//' @param diPixel  A vector of downstream pixels, diPixel[i] is downstream from pixel i
//' @param nx The total number of pixels in the network
//' @param value The value to be added (e.g., length) when computing distance
//' @return A matrix with dimensions `length(x)` by `nrow(ws)`
// [[Rcpp::export]]
NumericMatrix dmat(NumericVector x, NumericVector dsPixel, int nx, NumericVector value) {
	NumericMatrix distance(x.size(), nx);
	std::fill(distance.begin(), distance.end(), NumericVector::get_na());

	int count = 0;
	for(int xi = 0; xi < x.size(); ++xi) {
		int i = (int) x[xi] - 1; // -1 to correct for C++ indexing
		int j = i;
		double total = 0;
		distance(xi, j) = total;
		while(j > 0) {
			total += value(j);
			j = (int) dsPixel(j) - 1;
			distance(xi, j) = total;
			NumericVector::iterator inx = std::find(x.begin(), x.end(), j);
			if(inx != x.end())
				distance(inx - x.begin(), xi) = -1 * total;
			count++;
			if(count > nx*x.size())
				stop("Topology error");
		}
	}
	return distance;
}
