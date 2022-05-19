#include <Rcpp.h>
using namespace Rcpp;

//' Connect points in a watershed, accumulating a value as we go
//'
//' The values and ds pixels should match; i.e., values[i] is the value for pixel i, 
//' and dsPixel[i] is that pixel's downstream pixel
//'
//' @param dsPixel A vector of downstream pixels, diPixel[i] is downstream from pixel i
//' @param upstream The id of the upstream point
//' @param downstream The id of the downstream point
//' @param value A vector of values to accumulate
//' @return A list, the first entry is the connected pixels, the second the accumulation
//' up to (but excluding) each pixel
// [[Rcpp::export]]
List connectCPP(NumericVector dsPixel, int upstream, int downstream, 
		NumericVector value) {
	std::vector<int> connectedPts;
	std::vector<double> distance;
	double accum = 0.0;
	int current = upstream;
	connectedPts.push_back(current);
	distance.push_back(accum);
	while(current != downstream && !Rcpp::NumericVector::is_na(current)) {
		accum += value(current - 1);
		current = dsPixel(current - 1); // -1 to correct for C++ indexing
		connectedPts.push_back(current);
		distance.push_back(accum);
	}
	if(Rcpp::NumericVector::is_na(current) && !Rcpp::NumericVector::is_na(dsPixel(downstream))) {
		connectedPts.clear();
		distance.clear();
	}
	return Rcpp::List::create(connectedPts, distance);
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
