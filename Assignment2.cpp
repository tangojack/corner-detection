#include "stdAfx.h"
#include "Assignment2.h"
#include <cmath>
/////////////////////////////
// CCorner Source File //
/////////////////////////////
using namespace std;

// this function convert a given CImage to a GrayScale image
void CCorner::RGBToGrayScale(CImage* pIn, CImage* pOut)
{
	//
	// INPUT:
	//     CImage* pIn:		The input image with 24bit depth
	//
	// OUTPUT:
	//     CImage* pOut:	The output image. It has ALREADY been initialized
	//                      with the same dimension as the input image (pIN) and 
	//                      formatted to 8bit depth (256 gray levels). So, please
	//                      use 'SetIndex' instead of 'SetRGB'.
	//
	
	// Begin your code here: //
	byte* r = new byte();
	byte* g = new byte();
	byte* b = new byte();
	byte intensity;
	int height = pIn->GetHeight();
	int width = pIn->GetWidth();
	for (int i = 0; i < width; i++){
		for (int j = 0; j < height; j++){
			pIn->GetRGB(i, j, r, g, b);		
			intensity = (byte)((0.299 * (*r)) + (0.587 * (*g)) + (0.114 * (*b)));
			pOut->SetIndex(i, j, intensity);
		}
	}
}

// this function obtains the corners from a given GrayScale image
void CCorner::ObtainCorners(CImage* pIn, double sigma, double threshold, vector<C2DPoint*>* pResultCorners)
{
		//
	// INPUT:
	//     CImage* pIn:		The input grayscale image, which is exactly the output of
	//                      RGBToGrayScale function.
	//     double sigma:    The sigma value for your convolve function
	//     double threhold: The minimum value for qualifying a corner. Please refers to the lecture 3's notes
	//
	// OUTPUT:
	//     vector<C2DPoint*>* pResultCorners:	
	//                      A std::vector object that holds all detected corners. Each corner is represented by
	//                      a C2DPoint structure. An example is given below, showing how a corner object is
	//                      initialized and stored in pResultCorners:
	//                      
	//                      	C2DPoint* pPnt = new C2DPoint(x, y);
	//                      	pResultCorners.push_back(pPnt);
	//
	//
	// Begin your code here: //

	// Step 1: Compute a proper size for Gaussian filter
	
	int n = (int)ceil(sigma * sqrt(2*log(1000.0)));
	int size = 2 * n + 1;
	
	// Step 2: Define the Gaussian filter and partial filter
	double* mask = new double[size];	
	double sum = 0;
	int x;
	for(int i = 0; i < size; i++){
		x = i-n;
		mask[i]= exp(-(x*x)/(2*sigma*sigma));
		sum += mask[i];
	}
	for(int i = 0; i < size; i++){
		mask[i] /= sum;
	}

	
	// Step 3: Compute Ix, Iy 
	int height = pIn->GetHeight();
	int width = pIn->GetWidth();

	//Defining the double arrays that are going to be used and initialising each element to zero
	double **Ix = new double*[height];
	for (int i = 0; i < height; i++){
		Ix[i] = new double[width];
		for (int j = 0; j < width; j++){
			Ix[i][j] = 0;
		}
	}
	double **Iy = new double*[height];
	for (int i = 0; i < height; i++){
		Iy[i] = new double[width];
		for (int j = 0; j < width; j++){
			Iy[i][j] = 0;
		}
	}
	double **Ixx = new double*[height];
	for (int i = 0; i < height; i++){
		Ixx[i] = new double[width];
		for (int j = 0; j < width; j++){
			Ixx[i][j] = 0;
		}
	}
	double **Ixy = new double*[height];
	for (int i = 0; i < height; i++){
		Ixy[i] = new double[width];
		for (int j = 0; j < width; j++){
			Ixy[i][j] = 0;
		}
	}
	double **Iyy = new double*[height];
	for (int i = 0; i < height; i++){
		Iyy[i] = new double[width];
		for (int j = 0; j < width; j++){
			Iyy[i][j] = 0;
		}
	}

	//Smoothed array of Ixx
	double **smix = new double*[height];  
	for (int i = 0; i < height; i++){
		smix[i] = new double[width];
		for (int j = 0; j < width; j++){
			smix[i][j] = 0;
		}
	}

	//Smoothed array of Iyy
	double **smiy = new double*[height];
	for (int i = 0; i < height; i++){
		smiy[i] = new double[width];
		for (int j = 0; j < width; j++){
			smiy[i][j] = 0;
		}
	}

	//Smoothed array of Ixy
	double **smixy = new double*[height];
	for (int i = 0; i < height; i++){
		smixy[i] = new double[width];
		for (int j = 0; j < width; j++){
			smixy[i][j] = 0;
		}
	}

	//Smoothed arrays after only the Gaussian Filter in the x Direction
	double **xsmix = new double*[height];
	for (int i = 0; i < height; i++){
		xsmix[i] = new double[width];
		for (int j = 0; j < width; j++){
			xsmix[i][j] = 0;
		}
	}
	double **xsmiy = new double*[height];
	for (int i = 0; i < height; i++){
		xsmiy[i] = new double[width];
		for (int j = 0; j < width; j++){
			xsmiy[i][j] = 0;
		}
	}
	double **xsmixy = new double*[height];
	for (int i = 0; i < height; i++){
		xsmixy[i] = new double[width];
		for (int j = 0; j < width; j++){
			xsmixy[i][j] = 0;
		}
	}
	//Stores the R-Value
	double **r = new double*[height];
	for (int i = 0; i < height; i++){
		r[i] = new double[width];
		for (int j = 0; j < width; j++){
			r[i][j] = 0;
		}
	}
	int **maxima = new int*[height];
	for (int i = 0; i < height; i++){
		maxima[i] = new int[width];
		for (int j = 0; j < width; j++){
			maxima[i][j] = 0;
		}
	}


	//Computing Ix
	for (int y = 0; y < height; y++){ 
		Ix[y][0] = pIn->GetIndex(1, y) - pIn->GetIndex(0, y);
		for (int x = 1; x < width - 1; x++){
			Ix[y][x] = -(pIn->GetIndex(x-1, y) *0.5) + (pIn->GetIndex(x+1, y) * 0.5);
		}
		Ix[y][width - 1] = pIn->GetIndex(width - 1, y) - pIn->GetIndex(width - 2, y);
	}
	
	//Computing Iy
	for (int x = 0; x <width; x++){
		Iy[0][x] = pIn->GetIndex(x,1) - pIn->GetIndex(x,0);
		for (int y = 1; y < height - 1; y++){
			Iy[y][x] = (0.5 * pIn->GetIndex(x, y+1)) - (0.5 * pIn->GetIndex(x, y-1));
		}
		Iy[height - 1][x] = pIn->GetIndex(x, height - 1) - pIn->GetIndex(x, height - 2);
	}


	// Step 4: Compute Ixx, Iyy, Ixy
	
	for (int y = 0; y < height; y++){
		for (int x = 0; x < width; x++){
			Ixx[y][x] = Ix[y][x] * Ix[y][x];
			Iyy[y][x] = Iy[y][x] * Iy[y][x];
			Ixy[y][x] = Ix[y][x] * Iy[y][x];
		}
	}
	
	
	// Step 5: Smooth Ixx, Iyy, Ixy

	//Smoothing in Ixx in x direction
	double *temp;
	for (int y = 0; y < height; y++){
		for (int x = 0; x < width; x++){
			sum = 0;
			for (int i = -n; i <= n; i++){
				if (x + i < 0 || x + i >= width){
					continue;
				}
				sum += mask[i+n];
			}
			temp = new double[size];
			for (int i = -n; i <= n; i++){
				if (x + i < 0 || x + i >= width)
					temp[i + n] = 0;
				else
					temp[i + n] = mask[i + n] / sum;
			}
			xsmix[y][x] = 0;
			for (int i = -n; i <= n; i++){
				if (temp[i + n] == 0){
					continue;
				}	
				xsmix[y][x] += Ixx[y][x+i] * temp[i+n];
			}
		}
	}
	//Smoothing in y direction
	for (int x = 0; x < width; x++){
		for (int y = 0; y < height; y++){
			sum = 0;
			for (int i = -n; i <= n; i++){
				if (y + i < 0 || y + i >= height){
					continue;
				}
				sum += mask[i+n];
			}
			temp = new double[size];
			for (int i = -n; i <= n; i++){
				if (y + i < 0 || y + i >= height)
					temp[i + n] = 0;
				else
					temp[i + n] = mask[i + n] / sum;
			}
			smix[y][x] = 0;
			for (int i = -n; i <= n; i++){
				if (temp[i + n] == 0){
					continue;
				}	
				smix[y][x] += xsmix[y+i][x] * temp[i+n];
			}
		}
	}

	//Similarly, smoothing Iyy
	for (int y = 0; y < height; y++){
		for (int x = 0; x < width; x++){
			sum = 0;
			for (int i = -n; i <= n; i++){
				if (x + i < 0 || x + i >= width){
					continue;
				}
				sum += mask[i+n];
			}
			temp = new double[size];
			for (int i = -n; i <= n; i++){
				if (x + i < 0 || x + i >= width)
					temp[i + n] = 0;
				else
					temp[i + n] = mask[i + n] / sum;
			}
			xsmiy[y][x] = 0;
			for (int i = -n; i <= n; i++){
				if (temp[i + n] == 0){
					continue;
				}	
				xsmiy[y][x] += Iyy[y][x+i] * temp[i+n];
			}
		}
	}

	for (int x = 0; x < width; x++){
		for (int y = 0; y < height; y++){
			sum = 0;
			for (int i = -n; i <= n; i++){
				if (y + i < 0 || y + i >= height){
					continue;
				}
				sum += mask[i+n];
			}
			temp = new double[size];
			for (int i = -n; i <= n; i++){
				if (y + i < 0 || y + i >= height)
					temp[i + n] = 0;
				else
					temp[i + n] = mask[i + n] / sum;
			}
			smiy[y][x] = 0;
			for (int i = -n; i <= n; i++){
				if (temp[i + n] == 0){
					continue;
				}	
				smiy[y][x] += xsmiy[y+i][x] * temp[i+n];
			}
		}
	}

	//Smoothing Ixy
	for (int y = 0; y < height; y++){
		for (int x = 0; x < width; x++){
			sum = 0;
			for (int i = -n; i <= n; i++){
				if (x + i < 0 || x + i >= width){
					continue;
				}
				sum += mask[i+n];
			}
			temp = new double[size];
			for (int i = -n; i <= n; i++){
				if (x + i < 0 || x + i >= width)
					temp[i + n] = 0;
				else
					temp[i + n] = mask[i + n] / sum;
			}
			xsmixy[y][x] = 0;
			for (int i = -n; i <= n; i++){
				if (temp[i + n] == 0){
					continue;
				}	
				xsmixy[y][x] += Ixy[y][x+i] * temp[i+n];
			}
		}
	}

	for (int x = 0; x < width; x++){
		for (int y = 0; y < height; y++){
			sum = 0;
			for (int i = -n; i <= n; i++){
				if (y + i < 0 || y + i >= height){
					continue;
				}
				sum += mask[i+n];
			}
			temp = new double[size];
			for (int i = -n; i <= n; i++){
				if (y + i < 0 || y + i >= height)
					temp[i + n] = 0;
				else
					temp[i + n] = mask[i + n] / sum;
			}
			smixy[y][x] = 0;
			for (int i = -n; i <= n; i++){
				if (temp[i + n] == 0){
					continue;
				}	
				smixy[y][x] += xsmixy[y+i][x] * temp[i+n];
			}
		}
	}

	
	// Step 6: Compute R	
	for (int x = 0; x < width; x++){
		for (int y = 0; y < height; y++){
			r[y][x]	= (smix[y][x] * smiy[y][x]) - (smixy[y][x] * smixy[y][x]) - (0.04 * (smix[y][x] + smiy[y][x])); 	
		}
	}
	
	// Step 7: Locate Maxima in R
	int flag = 1;
	for (int x = 1; x < width - 1; x++){
		for (int y = 1; y < height - 1; y++){
			flag = 1;
			for (int i = -1; i <= 1; i++){
				for (int j = -1; j <= 1; j++){
					if (r[y+i][x+j] >= r[y][x] && !(i == 0 && j ==0)){  //Comparing neighbours
						flag = 0; 
					}
				}
			}
			if (flag == 1){
				maxima[y][x] = 1;
			}
		}
	}

	// Step 8: Compute corner candidates up to sub-pixel accuracy and interpolate R value for corner candidates.
	double subx, suby;
	double ax, bx, cx, ay, by, cy, rx, ry, rmax;
	for (int y = 0; y < height; y++){
		for (int x = 0; x < width; x++){
			if(maxima[y][x] == 1){
				ax = (r[y][x+1] + r[y][x-1]-(2.0*r[y][x])) / 2.0;
				bx = (r[y][x+1] - r[y][x-1]) / 2.0;
				cx = r[y][x];
				ay = (r[y+1][x] + r[y-1][x]-(2.0*r[y][x])) / 2.0;
				by = (r[y+1][x] - r[y-1][x]) / 2.0;
				cy = r[y][x]; 
				
				subx = -bx / (2.0 * ax);
				suby = -by / (2.0 * ay);
				rx = (ax * subx * subx) + (bx * subx) + cx;
				ry = (ay * suby * suby) + (by * suby) + cy;
				rmax = rx > ry ? rx : ry;
				if (rmax > threshold && rmax == ry){
					C2DPoint* pPnt = new C2DPoint(x, height - y - 1 + suby);
					(*pResultCorners).push_back(pPnt);
				}
				else if (rmax > threshold && rmax == rx){
					C2DPoint* pPnt = new C2DPoint(x+subx, height - y - 1);
					(*pResultCorners).push_back(pPnt);
				}
			}
		}
	}
	//Step 9 in done in the above step itself.
	// Step 9: Use the threshold value to identify strong corners for output
	

	//Freeing up memory
	for(int i = 0; i < height; ++i) {
        delete[] Ix[i];   
    }
    delete[] Ix;

    for(int i = 0; i < height; ++i) {
        delete[] Iy[i];   
    }
    delete[] Iy;

	for(int i = 0; i < height; ++i) {
        delete[] Ixx[i];   
    }
    delete[] Ixx;

    for(int i = 0; i < height; ++i) {
        delete[] Iyy[i];   
    }
    delete[] Iyy;
    
    for(int i = 0; i < height; ++i) {
        delete[] Ixy[i];   
    }
    delete[] Ixy;
    
    for(int i = 0; i < height; ++i) {
        delete[] smix[i];   
    }
    delete[] smix;
    
    for(int i = 0; i < height; ++i) {
        delete[] smixy[i];   
    }
    delete[] smixy;
    
    for(int i = 0; i < height; ++i) {
        delete[] smiy[i];   
    }
    delete[] smiy;
    
    for(int i = 0; i < height; ++i) {
        delete[] xsmix[i];   
    }
    delete[] xsmix;

    for(int i = 0; i < height; ++i) {
        delete[] xsmiy[i];   
    }
    delete[] xsmiy;

    for(int i = 0; i < height; ++i) {
        delete[] xsmixy[i];   
    }
    delete[] xsmixy;

    for(int i = 0; i < height; ++i) {
        delete[] r[i];   
    }
    delete[] r;

    for(int i = 0; i < height; ++i) {
        delete[] maxima[i];   
    }
    delete[] maxima;
}