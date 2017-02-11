//
// detect.cpp : Detect cars in satellite images.
//
// Based on skeleton code by D. Crandall, Spring 2017
//
// PUT YOUR NAMES HERE
//
//

#include <SImage.h>
#include <SImageIO.h>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

// The simple image class is called SDoublePlane, with each pixel represented as
// a double (floating point) type. This means that an SDoublePlane can represent
// values outside the range 0-255, and thus can represent squared gradient magnitudes,
// harris corner scores, etc. 
//
// The SImageIO class supports reading and writing PNG files. It will read in
// a color PNG file, convert it to grayscale, and then return it to you in 
// an SDoublePlane. The values in this SDoublePlane will be in the range [0,255].
//
// To write out an image, call write_png_file(). It takes three separate planes,
// one for each primary color (red, green, blue). To write a grayscale image,
// just pass the same SDoublePlane for all 3 planes. In order to get sensible
// results, the values in the SDoublePlane should be in the range [0,255].
//

// Below is a helper functions that overlays rectangles
// on an image plane for visualization purpose. 

// Draws a rectangle on an image plane, using the specified gray level value and line width.
//
void overlay_rectangle(SDoublePlane &input, int _top, int _left, int _bottom, int _right, double graylevel, int width)
{
  for(int w=-width/2; w<=width/2; w++) {
    int top = _top+w, left = _left+w, right=_right+w, bottom=_bottom+w;

    // if any of the coordinates are out-of-bounds, truncate them 
    top = min( max( top, 0 ), input.rows()-1);
    bottom = min( max( bottom, 0 ), input.rows()-1);
    left = min( max( left, 0 ), input.cols()-1);
    right = min( max( right, 0 ), input.cols()-1);
      
    // draw top and bottom lines
    for(int j=left; j<=right; j++)
	  input[top][j] = input[bottom][j] = graylevel;
    // draw left and right lines
    for(int i=top; i<=bottom; i++)
	  input[i][left] = input[i][right] = graylevel;
  }
}

// DetectedBox class may be helpful!
//  Feel free to modify.
//
class DetectedBox {
public:
  int row, col, width, height;
  double confidence;
};

// Function that outputs the ascii detection output file
void  write_detection_txt(const string &filename, const vector<DetectedBox> &cars)
{
  ofstream ofs(filename.c_str());

  for(vector<DetectedBox>::const_iterator s=cars.begin(); s != cars.end(); ++s)
    ofs << s->row << " " << s->col << " " << s->width << " " << s->height << " " << s->confidence << endl;
}

// Function that outputs a visualization of detected boxes
void  write_detection_image(const string &filename, const vector<DetectedBox> &cars, const SDoublePlane &input)
{
  SDoublePlane output_planes[3];

  for(int p=0; p<3; p++)
    {
      output_planes[p] = input;
      for(vector<DetectedBox>::const_iterator s=cars.begin(); s != cars.end(); ++s)
	overlay_rectangle(output_planes[p], s->row, s->col, s->row+s->height-1, s->col+s->width-1, p==2?255:0, 2);
    }

  SImageIO::write_png_file(filename.c_str(), output_planes[0], output_planes[1], output_planes[2]);
}



// The rest of these functions are incomplete. These are just suggestions to 
// get you started -- feel free to add extra functions, change function
// parameters, etc.

// Convolve an image with a separable convolution kernel
//
SDoublePlane convolve_separable(const SDoublePlane &input, const SDoublePlane &row_filter, const SDoublePlane &col_filter)
{
  SDoublePlane output_1D(input.rows(), input.cols());
  SDoublePlane output_final(input.rows(), input.cols());

  SDoublePlane input_convolve(input.rows()+2, input.cols()+2);
  int HEIGHT = input.rows();
  int WIDTH = input.cols();
  int height_inputC = input_convolve.rows();
  int width_inputC = input_convolve.cols();
  
  for(int r = 1; r < HEIGHT-1 ; r++)
  {
   for(int c = 1; c < WIDTH-1 ; c++)
   {
      input_convolve[r][c] = input[r-1][c-1];
   }
  }
  
  for(int c=0; c < WIDTH; c++)
  {
     input_convolve[0][c+1] = input[0][c];
     input_convolve[height_inputC-1][c+1] = input[HEIGHT-1][c];
  }

  for(int r = 0; r < HEIGHT; r++)
  {
     input_convolve[r+1][0] = input[r][0];
     input_convolve[r+1][width_inputC-1] = input[r][WIDTH-1];
  }
  // Convolution code here
  
  for(int row = 1; row < height_inputC-1; row++)
  { 
    for(int col = 1; col < width_inputC-1; col++ ) 
    {
        float accumulation = 0;
        float weightsum = 0; 
        for(int i = -1; i <= 1; i++ )
        {
                    double k = input_convolve[row+i][col];
                    accumulation += k * row_filter[1-i][0];
                    weightsum += row_filter[1-i][0];
        }
        output_1D[row-1][col-1] = (double)(accumulation/weightsum);
    }
  }
  
                
  for(int r = 1; r < HEIGHT-1 ; r++)
  {
   for(int c = 1; c < WIDTH-1 ; c++)
   {
      input_convolve[r][c] = output_1D[r-1][c-1];
   }
  }
  
  for(int c=0; c < WIDTH; c++)
  {
     input_convolve[0][c+1] = output_1D[0][c];
     input_convolve[height_inputC-1][c+1] = output_1D[HEIGHT-1][c];
  }

  for(int r = 0; r < HEIGHT; r++)
  {
     input_convolve[r+1][0] = output_1D[r][0];
     input_convolve[r+1][width_inputC-1] = output_1D[r][WIDTH-1];
  }
  // Convolution code here
  
  for(int row = 1; row < height_inputC-1; row++)
  { 
    for(int col = 1; col < width_inputC-1; col++ ) 
    {
        float accumulation = 0;
        float weightsum = 0; 
        for(int i = -1; i <= 1; i++ )
        {
                    double k = input_convolve[row][col+i];
                    accumulation += k * col_filter[0][1-i];
                    weightsum += col_filter[0][1-i];
        }
        output_final[row-1][col-1] = (double)(accumulation/weightsum);
    }
  }
  return output_final;
}

// Convolve an image with a  convolution kernel
//
SDoublePlane convolve_general(const SDoublePlane &input, const SDoublePlane &filter)
{

  // Convolution code here
  SDoublePlane output(input.rows(), input.cols());
  SDoublePlane input_convolve(input.rows()+2, input.cols()+2);
  int HEIGHT = input.rows();
  int WIDTH = input.cols();
  int height_inputC = input_convolve.rows();
  int width_inputC = input_convolve.cols();
  
  for(int r = 1; r < HEIGHT-1 ; r++)
  {
   for(int c = 1; c < WIDTH-1 ; c++)
   {
      input_convolve[r][c] = input[r-1][c-1];
   }
  }
  
  for(int c=0; c < WIDTH; c++)
  {
     input_convolve[0][c+1] = input[0][c];
     input_convolve[height_inputC-1][c+1] = input[HEIGHT-1][c];
  }

  for(int r = 0; r < HEIGHT; r++)
  {
     input_convolve[r+1][0] = input[r][0];
     input_convolve[r+1][width_inputC-1] = input[r][WIDTH-1];
  }
  //Boundary pixel values:
  input_convolve[0][0] = input[0][0];
  input_convolve[0][width_inputC-1] = input[0][WIDTH-1];
  input_convolve[height_inputC-1][0] = input[HEIGHT-1][0];
  input_convolve[height_inputC-1][width_inputC-1] = input[HEIGHT-1][WIDTH-1];
  for(int row = 1; row < height_inputC-1; row++)
  { 
    for(int col = 1; col < width_inputC-1; col++ ) 
    {
        float accumulation = 0;
        float weightsum = 0; 
        for(int i = -1; i <= 1; i++ )
        {
                for(int j = -1; j <= 1; j++ )
            {
                    double k = input_convolve[row+i][col+j];
                    accumulation += k * filter[1-i][1-j];
                    weightsum += filter[1-i][1-j];
                }
        }
        output[row-1][col-1] = (double)(accumulation/weightsum);
    }
  }
  return output;
}
SDoublePlane reflectImage(const SDoublePlane &input)
{
  
  SDoublePlane input_convolve(input.rows()+2, input.cols()+2);
  int HEIGHT = input.rows();
  int WIDTH = input.cols();
  int height_inputC = input_convolve.rows();
  int width_inputC = input_convolve.cols();

  for(int r = 1; r < HEIGHT-1 ; r++)
  {
   for(int c = 1; c < WIDTH-1 ; c++)
   {
      input_convolve[r][c] = input[r-1][c-1];
   }
  }
  
  for(int c=0; c < WIDTH; c++)
  {
     input_convolve[0][c+1] = input[0][c];
     input_convolve[height_inputC-1][c+1] = input[HEIGHT-1][c];
  }

  for(int r = 0; r < HEIGHT; r++)
  {
     input_convolve[r+1][0] = input[r][0];
     input_convolve[r+1][width_inputC-1] = input[r][WIDTH-1];
  }
  //Boundary pixel values:
  input_convolve[0][0] = input[0][0];
  input_convolve[0][width_inputC-1] = input[0][WIDTH-1];
  input_convolve[height_inputC-1][0] = input[HEIGHT-1][0];
  input_convolve[height_inputC-1][width_inputC-1] = input[HEIGHT-1][WIDTH-1];

  //cout<<"new image created"<<endl;
  return input_convolve;
}

SDoublePlane outputGradientOperator(const SDoublePlane &input_convolve, const SDoublePlane &filter)
{
 int height_inputC = input_convolve.rows();
 int width_inputC = input_convolve.cols();
 SDoublePlane output = SDoublePlane(height_inputC-2, width_inputC-2);
  for(int row = 1; row < height_inputC-1; row++)
  { 
    for(int col = 1; col < width_inputC-1; col++ ) 
    {
        float accumulation = 0;
        float weightsum = 0; 
        for(int i = -1; i <= 1; i++ )
        {
                for(int j = -1; j <= 1; j++ )
            {
                    double k = input_convolve[row+i][col+j];
                    accumulation += k * filter[1+i][1+j];
                    //weightsum += filter[1+i][1+j];
                }
        }
        output[row-1][col-1] = (double)(accumulation);
    }
  }
  //cout<<"convolution done"<<endl;  
  return output;
}
// Found edge 
bool foundEdge(const SDoublePlane &output_after_nms, int window, int row, int col)
{
for (int r = row-window; r<=row+window; r++)
{
	for(int c = col-window; col<=col+window; c++)
	{
		if (output_after_nms[r][c] >= 150.0)
			return true;
	}
}
return false;

}

// Apply Hystersis thresholding

void hysterisisThresholding(const SDoublePlane &output_after_nms, const SDoublePlane &output_hysterisis)
{
int height = output_after_nms.rows();
int width = output_after_nms.cols();
double low = 70.0, high = 150.0;
 for (int r = 2; r< height - 2; r++)
 {
 	for (int c = 2; c<width-2; c++)
 	{
if (output_after_nms[r][c] >= high )
	//output_hysterisis[r][c] = output_after_nms[r][c];
		output_hysterisis[r][c] = 255.0;

else if (output_after_nms[r][c] > low && output_after_nms[r][c] <high)
{
bool found = false;

if (foundEdge(output_after_nms, 1, r, c))
{
	found = true;
	//output_hysterisis[r][c] = output_after_nms[r][c];
	output_hysterisis[r][c] = 255.0;
}
else  if (!found && foundEdge(output_after_nms, 2, r, c))
{
//output_hysterisis[r][c] = output_after_nms[r][c];
		output_hysterisis[r][c] = 255.0;

}


}

}
}
}

// Apply Non Max Suppression
void nonMaxSuppression(const SDoublePlane &output_after_nms, const SDoublePlane &output, const SDoublePlane &orientation)
{
 int height = output.rows();
 int width = output.cols();
 
 for (int r = 1; r< height - 1; r++)
 {
 	for (int c = 1; c<width-1; c++)
 	{

 	double orient = orientation[r][c];
 	double magnitude = output[r][c];
 	if ((orient == 0) && magnitude >= max(output[r][c-1], output[r][c+1]))
 		output_after_nms[r][c] = magnitude;
 	else if ((orient == 45) && magnitude >= max( output[r+1][c-1], output[r-1][c+1]))
 		output_after_nms[r][c] = magnitude;

 	else if ((orient == 90) && magnitude >= max (output[r-1][c], output[r+1][c] ))
 		output_after_nms[r][c] = magnitude;
    else if ((orient == 135) && magnitude >= max (output[r+1][c+1], output[r-1][c-1] ))
 		output_after_nms[r][c] = magnitude;
 	//cout << "output"<< output_after_nms[r][c];
 	//cout <<"\n";
 	
 	}
 }
}

// Apply a sobel operator to an image, returns the result
// 

SDoublePlane sobel_gradient_filter(const SDoublePlane &input, bool _gx)
{
  int rows = input.rows();
  int cols = input.cols();
  SDoublePlane output(rows, cols);
  SDoublePlane orientation(rows, cols);
  SDoublePlane output_after_nms = SDoublePlane(rows, cols);
  SDoublePlane output_hysterisis = SDoublePlane(rows, cols);
  // Implement a sobel gradient estimation filter with 1-d filters
  	SDoublePlane sobel_Gx = SDoublePlane(3,3);
	SDoublePlane sobel_Gy = SDoublePlane(3,3);
    sobel_Gx[0][0] = -1;
    sobel_Gx[0][1] = 0;
	sobel_Gx[0][2] = 1;
	sobel_Gx[1][0] = -2;
	sobel_Gx[1][1] = 0;
	sobel_Gx[1][2] = 2;
	sobel_Gx[2][0] = -1;
	sobel_Gx[2][1] = 0;
	sobel_Gx[2][2] = 1;
	sobel_Gy[0][0] = -1;
	sobel_Gy[0][1] = -2;
	sobel_Gy[0][2] = -1;
	sobel_Gy[1][0] = 0;
	sobel_Gy[1][1] = 0;
	sobel_Gy[1][2] = 0;
	sobel_Gy[2][0] = 1;
	sobel_Gy[2][1] = 2;
	sobel_Gy[2][2] = 1;
   
  SDoublePlane input_convolve = reflectImage(input);
  SDoublePlane output_Gx = outputGradientOperator(input_convolve, sobel_Gx);
 
  SDoublePlane output_Gy = outputGradientOperator(input_convolve, sobel_Gy);
  
  double pi = 180.0;
  for(int r = 1; r < rows-1; r++)
  {
   for(int c = 1; c < cols-1; c++)
   {
      output[r][c] = sqrt (pow(output_Gx[r][c] , 2) + pow(output_Gy[r][c] , 2));
      
      
     // output[r][c] = 0;  
      double tangX = output_Gx[r][c];
      double tangY = output_Gy[r][c];
      double temp = 0.0;
      temp = tangY/tangX;
      if (temp < 0)
      	temp = ((atan(temp)*180)/3.14 ) + 180;
      else
      	temp = ((atan(temp)*180)/3.14 );
      if (temp >22.5 && temp <=67.5)
      	temp = 45;
      else if (temp > 67.5 && temp <=112.5)
      	temp = 90;
       else if (temp > 112.5 && temp <=157.5)
      	temp = 135;
      else if ((temp >157.5 && temp <=202.5) || (temp > 0 && temp <=22.5) || (temp > 337.5 && temp <=360) )
      	temp = 0;
      else if (temp > 202.5 && temp <=247.5)
      	temp = 45;
      else if (temp >247.5 && temp <=292.5)
      	temp = 90;
      else
      	temp = 135;


      	orientation[r][c] = temp;
      	//cout<< orientation[r][c];
      	//cout <<"\n";
   }

  }

nonMaxSuppression(output_after_nms, output, orientation);

hysterisisThresholding(output_after_nms, output_hysterisis);


  return output_hysterisis;
}

// Apply an edge detector to an image, returns the binary edge map
// 
SDoublePlane find_edges(const SDoublePlane &input, double thresh=0)
{
   // Implement an edge detector of your choice, e.g.
  // use your sobel gradient operator to compute the gradient magnitude and threshold
  SDoublePlane output(input.rows(), input.cols());
  output = sobel_gradient_filter(input, true);
 

  return output;
}

bool isCar(int r, int c, SDoublePlane &outputEdges){
int countZeroes = 0;
int countWhites = 0;

//cout<< "in is car ";
//cout <<"\n";
//cout << r<<" "<<c<<endl;
for(int i = c; i< c+22; i++)
{
//count the number of zeroes(black pixels) on the first and last row
if (outputEdges[r][i] == 0)
countZeroes ++;
if (outputEdges[r+41][i] == 0)
countZeroes ++;
	
}
for(int i = r; i< r+42; i++)

{
//count the number of zeroes(black pixels) on the first and last row
if (outputEdges[i][c] == 0)
countZeroes ++;
if (outputEdges[i][c+21] == 0)
countZeroes ++;

}

if (countZeroes > 50)
{

for(int i = r+1; i< r+41; i++)
{
	for(int j = c+1; j< c+21; j++)
	{ 	
		if (outputEdges[i][j] == 255)
			countWhites ++;
		
	}
}
//cout << "Count of whites "<< countWhites<<endl;
if(countWhites > 85)
		return true;
}

return false;
}


//
// This main file just outputs a few test images. You'll want to change it to do 
//  something more interesting!
//
int main(int argc, char *argv[])
{
  if(!(argc == 2))
    {
      cerr << "usage: " << argv[0] << " input_image" << endl;
      return 1;
    }

  string input_filename(argv[1]);
  SDoublePlane input_image= SImageIO::read_png_file(input_filename.c_str());
  
  // test step 2 by applying mean filters to the input image
  SDoublePlane mean_filter(3,3);
  for(int i=0; i<3; i++)
    for(int j=0; j<3; j++)
      mean_filter[i][j] = 1/9.0;
  
	SDoublePlane gaussian_filter = SDoublePlane(3,3);
	gaussian_filter[0][0] = 0.0625;
	gaussian_filter[0][1] = 0.125;
	gaussian_filter[0][2] = 0.0625;
	gaussian_filter[1][0] = 0.125;
	gaussian_filter[1][1] = 0.25;
	gaussian_filter[1][2] = 0.125;
	gaussian_filter[2][0] = 0.0625;
	gaussian_filter[2][1] = 0.125;
	gaussian_filter[2][2] = 0.0625;
  SDoublePlane output_image = convolve_general(input_image, gaussian_filter);

  SImageIO::write_png_file("ConvolveGeneral.png",output_image,output_image,output_image);

  SDoublePlane separable_Hx = SDoublePlane(3,1);
  SDoublePlane separable_Hy = SDoublePlane(1,3);
  separable_Hx[0][0] = 0.25;
  separable_Hx[1][0] = 0.5;
  separable_Hx[2][0] = 0.25;

  separable_Hy[0][0] = 0.25;
  separable_Hy[0][1] = 0.5;
  separable_Hy[0][2] = 0.25;

  SDoublePlane output_image_sep = convolve_separable(input_image, separable_Hx, separable_Hy);
  
  SImageIO::write_png_file("ConvolveSeparable.png",output_image_sep,output_image_sep,output_image_sep);
 int count = 0;
 int height = input_image.rows();
 int width = input_image.cols();
    SDoublePlane output_image_demo = SDoublePlane(height, width);
 while(count < 10)
 { 
 output_image_demo = convolve_general(input_image, gaussian_filter);
    input_image = output_image_demo;
 	count ++;
 }
  SDoublePlane outputGxSobel = sobel_gradient_filter(output_image_demo, true);

  SImageIO::write_png_file("ConvolveSobelGx.png",outputGxSobel,outputGxSobel,outputGxSobel);

  // randomly generate some detected cars -- you'll want to replace this
  //  with your car detection code obviously!
bool found = false;
  SDoublePlane outputEdges = find_edges(input_image);
  vector<DetectedBox> cars;
  for(int r=0; r<height-42; r++)
    {
    found = false;
    for(int c=0; c<width-22; c++){
      DetectedBox s;
   	//cout <<"Check for the car"<<endl;
   	//cout << r<<" "<<c<<endl;
      if (isCar(r, c, outputEdges)){
      	s.row = r+1;
      	s.col = c+1;
     	s.width = 20;
      	s.height = 40;
      	//cout << "add to car"<<endl;
      	found = true;
      	c += 20;
      	s.confidence = rand();
      	cars.push_back(s);
    	 }

 	}
 	//cout<<" Found CARSSSSSSSSSSSSSSSSSSSSSSSSS"<<endl;
 	if (found)
 		r += 40;
    
    	}
  write_detection_txt("detected.txt", cars);
  write_detection_image("detected.png", cars, input_image);
}
